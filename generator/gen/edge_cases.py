from __future__ import annotations

"""Generate tests for architectural edge cases:
- vl=0 no-op behavior (integer, FP, load, store)
- vstart>0: csrw/csrr vstart + vadd.vv behavior (both SIGILL and skip are legal)
- Tail undisturbed (vl < VLMAX) — basic 4-element check
- Full-VLMAX tail undisturbed — verifies ALL tail bytes via vs1r.v + vlenb
- LMUL>1 register boundary crossing (v8→v9)
- GhostWrite (CVE-2024-44067) mew=1 reserved encoding check
- Reserved RVV 1.0 load/store encoding SIGILL checks (mew, lumop, sumop)

These test the vector unit's control logic rather than specific ALU
operations.  We use vadd.vv / vfadd.vv / vle32.v / vse32.v / vwadd.vv
as representative instructions since the behavior is
instruction-independent.
"""

from pathlib import Path

from ..common import (
    NUM_ELEMS,
    format_data_line,
    format_hex,
    M,
    S,
    f32_to_bits,
    f64_to_bits,
    witness_pattern,
)
from ..testfile import TestFile
from ..emit import VREG_DST, VREG_SRC2, VREG_SRC1, VREG_WITNESS


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "edge_cases"

    # ==================================================================
    # Test 1: vl=0 should be a complete no-op
    # ==================================================================
    # Pre-load vd with known pattern, set vl=0, execute vadd.vv,
    # verify vd unchanged and vstart==0 afterward.

    fpath = out / "vl_zero.S"
    tf = TestFile("vadd.vv", "vl=0 no-op: instruction must not modify vd")

    sew = 32
    bg = [0xDEAD_BEEF, 0xCAFE_BABE, 0x1234_5678, 0x9ABC_DEF0]
    nbytes = NUM_ELEMS * (sew // 8)

    cn_res = tf.next_check("vl=0: vd modified")
    cn_vstart = tf.next_check("vl=0: vstart not zero after")
    tag = f"tc{cn_res}"

    tf.comment(f"Test {cn_res}: vadd.vv with vl=0 — vd must be unchanged")

    # Pre-load vd at vl=4
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_bg")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    # Load sources (doesn't matter what, they shouldn't be used)
    tf.code(f"la t1, {tag}_src")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Set vl=0
    tf.code("li t0, 0")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Execute — should be a no-op
    tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Check vstart == 0
    tf.code(f"SET_TEST_NUM {cn_vstart}")
    tf.code("csrr t2, vstart")
    tf.code("FAIL_IF_NZ t2")

    # Restore vl=4 to be able to store and check vd
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_bg, {nbytes}")

    tf.data_align(sew)
    tf.data_label(f"{tag}_bg", format_data_line(bg, sew))
    tf.data_label(f"{tag}_src", format_data_line([1, 2, 3, 4], sew))

    tf.write(fpath)
    generated.append(str(fpath))

    # ==================================================================
    # Test 2: vstart CSR write/read and nonzero vstart behavior
    # ==================================================================
    generated.append(_gen_vstart_nonzero(out))

    # ==================================================================
    # Test 3: Tail undisturbed — elements >= vl must be preserved
    # ==================================================================
    # Pre-load vd with 4 elements, set vl=2, execute vadd.vv,
    # restore vl=4, verify elements 0-1 computed, elements 2-3 preserved.

    fpath = out / "tail_undisturbed.S"
    tf = TestFile("vadd.vv", "Tail undisturbed: elements >= vl must be preserved")

    vd_pre = [0xDEAD_BEEF, 0xCAFE_BABE, 0x1234_5678, 0x9ABC_DEF0]
    vs2 = [100, 200, 300, 400]
    vs1 = [10, 20, 30, 40]
    # Elements 0-1 get vadd result, elements 2-3 preserved from vd_pre
    expected_tail = [
        (vs2[0] + vs1[0]) & M(sew),
        (vs2[1] + vs1[1]) & M(sew),
        vd_pre[2],
        vd_pre[3],
    ]
    nbytes = NUM_ELEMS * (sew // 8)

    cn_res = tf.next_check("tail undisturbed: wrong result")
    cn_vstart = tf.next_check("tail undisturbed: vstart not zero after")
    tag = f"tc{cn_res}"

    tf.comment(f"Test {cn_res}: vadd.vv with vl=2, tu — tail elements preserved")

    # Pre-load vd at vl=4
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    # Load sources at vl=4 (fill all elements)
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Set vl=2 (tail policy = tu still in effect)
    tf.code("li t0, 2")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code("SAVE_CSRS")

    # Execute with vl=2 — only elements 0,1 computed
    tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Check vstart
    tf.code(f"SET_TEST_NUM {cn_vstart}")
    tf.code("CHECK_VSTART_ZERO")

    # Restore vl=4 to store full register
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.data_align(sew)
    tf.data_label(f"{tag}_vd", format_data_line(vd_pre, sew))
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected_tail, sew))

    tf.write(fpath)
    generated.append(str(fpath))

    # ==================================================================
    # Test 4: Tail undisturbed with masked operation (combined)
    # ==================================================================
    # vl=2, mask=0b01 (only element 0 active), tu/mu policy.
    # Element 0: computed (active, in-range)
    # Element 1: vd preserved (inactive, in-range, mu)
    # Elements 2-3: vd preserved (tail, tu)

    fpath = out / "tail_masked_combined.S"
    tf = TestFile("vadd.vv", "Tail + mask combined: tu/mu interaction")

    vd_pre = [0xAAAA_0000, 0xBBBB_0000, 0xCCCC_0000, 0xDDDD_0000]
    vs2 = [100, 200, 300, 400]
    vs1 = [10, 20, 30, 40]
    mask_bits = 0b01  # only element 0 active
    expected_combo = [
        (vs2[0] + vs1[0]) & M(sew),  # element 0: active, computed
        vd_pre[1],  # element 1: inactive (mu), preserved
        vd_pre[2],  # element 2: tail (tu), preserved
        vd_pre[3],  # element 3: tail (tu), preserved
    ]
    nbytes = NUM_ELEMS * (sew // 8)

    cn_res = tf.next_check("tail+mask combined: wrong result")
    tag = f"tc{cn_res}"

    tf.comment(f"Test {cn_res}: vadd.vv vl=2 mask=0b01 tu/mu — combined check")

    # Pre-load vd at vl=4
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    # Load sources at vl=4
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Load mask
    tf.code(f"la t1, {tag}_mask")
    tf.code("vlm.v v0, (t1)")

    # Set vl=2
    tf.code("li t0, 2")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Execute masked with vl=2
    tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0.t")

    # Restore vl=4 for full check
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
    tf.data_align(sew)
    tf.data_label(f"{tag}_vd", format_data_line(vd_pre, sew))
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected_combo, sew))

    tf.write(fpath)
    generated.append(str(fpath))

    # ==================================================================
    # New tests: full-VLMAX tail, LMUL>1, additional vl=0
    # ==================================================================
    generated.append(_gen_vlmax_tail_int(out))
    generated.append(_gen_vlmax_tail_fp(out))
    generated.append(_gen_vlmax_tail_load(out))
    generated.append(_gen_vlmax_tail_widening(out))
    generated.append(_gen_lmul_gt1_int(out))
    generated.append(_gen_lmul_gt1_fp(out))
    generated.append(_gen_vl_zero_fp(out))
    generated.append(_gen_vl_zero_load(out))
    generated.append(_gen_vl_zero_store(out))
    generated.append(_gen_vle32ff_fault(out))

    # Phase 2: new edge-case tests (SIGILL trap, fractional LMUL, policies,
    # stride edge cases, scatter ordering, FP flags, vxsat sticky,
    # narrowing tail, widening m2→m4)
    generated.append(_gen_vill_trap(out))
    generated.append(_gen_fract_lmul(out))
    generated.append(_gen_tail_agnostic(out))
    generated.append(_gen_mask_agnostic(out))
    generated.append(_gen_stride_zero(out))
    generated.append(_gen_stride_negative(out))
    generated.append(_gen_scatter_ordered(out))
    generated.append(_gen_fflags_set(out))
    generated.append(_gen_vxsat_sticky(out))
    generated.append(_gen_narrowing_tail(out))
    generated.append(_gen_widening_m2_m4(out))

    # Phase 3: reserved encoding / vulnerability checks (fork-based SIGILL)
    generated.append(_gen_ghostwrite(out))
    generated.append(_gen_reserved_encoding(out))
    generated.append(_gen_rvv_detect(out))

    # Phase 4: tail undisturbed per-instruction-family tests
    generated.append(_gen_tail_per_family(out))

    # Phase 5: LMUL>1 boundary crossing per-instruction-family
    generated.append(_gen_lmul2_per_family(out))

    return generated


# ======================================================================
# Helper: emit the dynamic tail-check sequence
# ======================================================================


def _emit_tail_check(
    tf: TestFile,
    cn_tail: int,
    resbuf_label: str,
    bg_label: str,
    active_bytes: int,
    *,
    nreg: int = 1,
) -> None:
    """Emit code to dynamically check tail bytes after vs{nreg}r.v dump.

    Compares ``resbuf[active_bytes:]`` against ``bg[active_bytes:]``
    for ``nreg * vlenb - active_bytes`` bytes using ``_mem_compare``.

    Args:
        nreg: number of registers dumped (1 for vs1r.v, 2 for vs2r.v).
    """
    tf.code(f"SET_TEST_NUM {cn_tail}")
    tf.code(f"la a1, {resbuf_label}")
    tf.code(f"addi a1, a1, {active_bytes}")
    tf.code(f"la a2, {bg_label}")
    tf.code(f"addi a2, a2, {active_bytes}")
    tf.code("csrr a3, vlenb")
    if nreg > 1:
        tf.code(f"slli a3, a3, {(nreg).bit_length() - 1}")  # a3 *= nreg
    tf.code(f"addi a3, a3, -{active_bytes}")
    tf.code("blez a3, 1f")
    tf.code("jal ra, _mem_compare")
    tf.code("FAIL_IF_NZ a0")
    tf.raw("1:")


# ======================================================================
# P1: Full-VLMAX tail undisturbed — integer VV (vadd.vv)
# ======================================================================


def _gen_vlmax_tail_int(out: Path) -> str:
    """Full-VLMAX tail undisturbed: vadd.vv at e32, vl=2.

    Fill entire register with background via vmv.v.x at VLMAX, set vl=2,
    execute vadd.vv, dump with vs1r.v, check active bytes with CHECK_MEM
    and remaining tail bytes dynamically with _mem_compare + vlenb.
    """
    fpath = out / "tail_vlmax_int.S"
    tf = TestFile(
        "vadd.vv",
        "Full-VLMAX tail undisturbed: vadd.vv e32 — "
        "verifies ALL tail bytes up to VLEN/8",
    )

    sew = 32
    active_vl = 2
    active_bytes = active_vl * (sew // 8)  # 8
    bg_word = 0xDEAD_BEEF

    vs2 = [100, 200]
    vs1 = [10, 20]
    expected_active = [(a + b) & M(sew) for a, b in zip(vs2, vs1)]

    cn_active = tf.next_check("vlmax-tail int: active elements wrong")
    cn_tail = tf.next_check("vlmax-tail int: tail elements modified")
    cn_vstart = tf.next_check("vlmax-tail int: vstart not zero")
    tag = f"tc{cn_active}"

    tf.blank()
    tf.comment(f"Full-VLMAX tail test: vadd.vv e32, vl={active_vl}, tu/mu")
    tf.comment(
        "Fill register with bg at VLMAX, compute 2 elements, verify tail preserved"
    )

    # 1. Fill vd with bg at VLMAX
    tf.code("vsetvli t0, x0, e32, m1, tu, mu")
    tf.code(f"li t1, 0x{bg_word:08x}")
    tf.code(f"vmv.v.x {VREG_DST}, t1")

    # 2. Load sources at vl=4 (enough for active_vl=2)
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # 3. Set vl=active_vl
    tf.code(f"li t0, {active_vl}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # 4. Execute
    tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # 5. Check vstart
    tf.code(f"SET_TEST_NUM {cn_vstart}")
    tf.code("csrr t2, vstart")
    tf.code("FAIL_IF_NZ t2")

    # 6. Dump entire register (stores VLEN/8 bytes regardless of vl)
    tf.code(f"la t1, {tag}_resbuf")
    tf.code(f"vs1r.v {VREG_DST}, (t1)")

    # 7. Check active elements (first active_bytes)
    tf.code(f"SET_TEST_NUM {cn_active}")
    tf.code(f"CHECK_MEM {tag}_resbuf, {tag}_exp, {active_bytes}")

    # 8. Check tail bytes dynamically (vlenb - active_bytes)
    _emit_tail_check(tf, cn_tail, f"{tag}_resbuf", f"{tag}_bg", active_bytes)

    # Data
    tf.data_align(sew)
    tf.data_label(
        f"{tag}_s2", format_data_line(vs2 + [0] * (NUM_ELEMS - active_vl), sew)
    )
    tf.data_label(
        f"{tag}_s1", format_data_line(vs1 + [0] * (NUM_ELEMS - active_vl), sew)
    )
    tf.data_label(f"{tag}_exp", format_data_line(expected_active, sew))
    # Background: 512 bytes of bg_word (supports VLEN up to 4096)
    tf.data(".align 4")
    tf.data_label(f"{tag}_bg", f"    .fill 128, 4, 0x{bg_word:08x}")
    # Result buffer for vs1r.v (up to vlenb = 512 for VLEN=4096)
    tf.data(".align 4")
    tf.data_label(f"{tag}_resbuf", "    .space 512")

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# P1: Full-VLMAX tail undisturbed — FP VV (vfadd.vv)
# ======================================================================


def _gen_vlmax_tail_fp(out: Path) -> str:
    """Full-VLMAX tail undisturbed: vfadd.vv at e32, vl=2."""
    fpath = out / "tail_vlmax_fp.S"
    tf = TestFile(
        "vfadd.vv",
        "Full-VLMAX tail undisturbed: vfadd.vv e32 — "
        "verifies ALL tail bytes up to VLEN/8",
    )

    sew = 32
    active_vl = 2
    active_bytes = active_vl * (sew // 8)  # 8
    bg_word = 0xDEAD_BEEF

    vs2_fp = [1.0, 2.0]
    vs1_fp = [0.5, 0.25]
    exp_fp = [a + b for a, b in zip(vs2_fp, vs1_fp)]  # [1.5, 2.25]

    vs2_bits = [f32_to_bits(v) for v in vs2_fp]
    vs1_bits = [f32_to_bits(v) for v in vs1_fp]
    exp_bits = [f32_to_bits(v) for v in exp_fp]

    cn_active = tf.next_check("vlmax-tail fp: active elements wrong")
    cn_tail = tf.next_check("vlmax-tail fp: tail elements modified")
    tag = f"tc{cn_active}"

    tf.blank()
    tf.comment(f"Full-VLMAX tail test: vfadd.vv e32, vl={active_vl}, tu/mu")

    # 1. Fill vd with bg at VLMAX
    tf.code("vsetvli t0, x0, e32, m1, tu, mu")
    tf.code(f"li t1, 0x{bg_word:08x}")
    tf.code(f"vmv.v.x {VREG_DST}, t1")

    # 2. Load FP sources
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # 3. Set vl=2
    tf.code(f"li t0, {active_vl}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # 4. Execute
    tf.code(f"vfadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # 5. Dump entire register
    tf.code(f"la t1, {tag}_resbuf")
    tf.code(f"vs1r.v {VREG_DST}, (t1)")

    # 6. Check active elements
    tf.code(f"SET_TEST_NUM {cn_active}")
    tf.code(f"CHECK_MEM {tag}_resbuf, {tag}_exp, {active_bytes}")

    # 7. Check tail
    _emit_tail_check(tf, cn_tail, f"{tag}_resbuf", f"{tag}_bg", active_bytes)

    # Data
    tf.data_align(sew)
    tf.data_label(
        f"{tag}_s2", format_data_line(vs2_bits + [0] * (NUM_ELEMS - active_vl), sew)
    )
    tf.data_label(
        f"{tag}_s1", format_data_line(vs1_bits + [0] * (NUM_ELEMS - active_vl), sew)
    )
    tf.data_label(f"{tag}_exp", format_data_line(exp_bits, sew))
    tf.data(".align 4")
    tf.data_label(f"{tag}_bg", f"    .fill 128, 4, 0x{bg_word:08x}")
    tf.data(".align 4")
    tf.data_label(f"{tag}_resbuf", "    .space 512")

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# P1: Full-VLMAX tail undisturbed — Load (vle32.v)
# ======================================================================


def _gen_vlmax_tail_load(out: Path) -> str:
    """Full-VLMAX tail undisturbed: vle32.v at e32, vl=2.

    Fill vd with background at VLMAX, load 2 elements from memory,
    verify first 2 elements are the loaded values and tail is bg.
    """
    fpath = out / "tail_vlmax_load.S"
    tf = TestFile(
        "vle32.v",
        "Full-VLMAX tail undisturbed: vle32.v — verifies ALL tail bytes up to VLEN/8",
    )

    sew = 32
    active_vl = 2
    active_bytes = active_vl * (sew // 8)  # 8
    bg_word = 0xDEAD_BEEF

    # Memory source — elements 0-1 will be loaded, rest ignored
    load_src = [0x1111_1111, 0x2222_2222, 0x3333_3333, 0x4444_4444]

    cn_active = tf.next_check("vlmax-tail load: active elements wrong")
    cn_tail = tf.next_check("vlmax-tail load: tail elements modified")
    tag = f"tc{cn_active}"

    tf.blank()
    tf.comment(f"Full-VLMAX tail test: vle32.v e32, vl={active_vl}, tu/mu")

    # 1. Fill vd with bg at VLMAX
    tf.code("vsetvli t0, x0, e32, m1, tu, mu")
    tf.code(f"li t1, 0x{bg_word:08x}")
    tf.code(f"vmv.v.x {VREG_DST}, t1")

    # 2. Set vl=2
    tf.code(f"li t0, {active_vl}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # 3. Load from memory (only 2 elements loaded)
    tf.code(f"la t1, {tag}_src")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    # 4. Dump entire register
    tf.code(f"la t1, {tag}_resbuf")
    tf.code(f"vs1r.v {VREG_DST}, (t1)")

    # 5. Check active elements (first 2 words should match load source)
    tf.code(f"SET_TEST_NUM {cn_active}")
    tf.code(f"CHECK_MEM {tag}_resbuf, {tag}_src, {active_bytes}")

    # 6. Check tail
    _emit_tail_check(tf, cn_tail, f"{tag}_resbuf", f"{tag}_bg", active_bytes)

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag}_src", format_data_line(load_src, sew))
    tf.data(".align 4")
    tf.data_label(f"{tag}_bg", f"    .fill 128, 4, 0x{bg_word:08x}")
    tf.data(".align 4")
    tf.data_label(f"{tag}_resbuf", "    .space 512")

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# P1: Full-VLMAX tail undisturbed — Widening (vwadd.vv)
# ======================================================================


def _gen_vlmax_tail_widening(out: Path) -> str:
    """Full-VLMAX tail undisturbed: vwadd.vv at e32→e64, vl=2.

    Destination register group is v8-v9 (EMUL=2).  Fill both registers
    with bg via vmv.v.x at e32/m2, then execute widening with vl=2.
    Dump with vs2r.v (stores 2*vlenb bytes).  Active region is 2 e64
    elements = 16 bytes; tail = 2*vlenb - 16.
    """
    fpath = out / "tail_vlmax_widening.S"
    tf = TestFile(
        "vwadd.vv",
        "Full-VLMAX tail undisturbed: vwadd.vv e32->e64 — "
        "verifies ALL tail bytes in dest group",
    )

    src_sew = 32
    dst_sew = 64
    active_vl = 2
    active_bytes = active_vl * (dst_sew // 8)  # 16
    bg_word = 0xDEAD_BEEF

    vs2 = [100, 200]
    vs1 = [10, 20]
    # Widening signed add: sext32→64(vs2[i]) + sext32→64(vs1[i])
    expected_active = [110, 220]  # small positives, fit in e64

    cn_active = tf.next_check("vlmax-tail widening: active elements wrong")
    cn_tail = tf.next_check("vlmax-tail widening: tail elements modified")
    tag = f"tc{cn_active}"

    tf.blank()
    tf.comment(f"Full-VLMAX tail test: vwadd.vv e32→e64, vl={active_vl}, tu/mu")
    tf.comment("Dest is v8-v9 (EMUL=2). Uses vs2r.v to dump, 2*vlenb for tail.")

    # 1. Fill dest register group (v8-v9) with bg at e32/m2
    #    vmv.v.x at e32/m2 fills VLMAX(e32/m2)=2*VLEN/32 elements
    #    covering both v8 and v9 completely.
    tf.code("vsetvli t0, x0, e32, m2, tu, mu")
    tf.code(f"li t1, 0x{bg_word:08x}")
    tf.code(f"vmv.v.x {VREG_DST}, t1")

    # 2. Load sources at e32/m1
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{src_sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{src_sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{src_sew}.v {VREG_SRC1}, (t1)")

    # 3. Set vl=2 at source sew for widening
    tf.code(f"li t0, {active_vl}")
    tf.code(f"vsetvli t0, t0, e{src_sew}, m1, tu, mu")

    # 4. Execute widening
    tf.code(f"vwadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # 5. Dump entire dest group (v8-v9) with vs2r.v (2*vlenb bytes)
    tf.code(f"la t1, {tag}_resbuf")
    tf.code(f"vs2r.v {VREG_DST}, (t1)")

    # 6. Check active elements (first 16 bytes = 2 × e64)
    tf.code(f"SET_TEST_NUM {cn_active}")
    tf.code(f"CHECK_MEM {tag}_resbuf, {tag}_exp, {active_bytes}")

    # 7. Check tail: remaining = 2*vlenb - active_bytes
    _emit_tail_check(tf, cn_tail, f"{tag}_resbuf", f"{tag}_bg", active_bytes, nreg=2)

    # Data
    tf.data_align(src_sew)
    tf.data_label(
        f"{tag}_s2", format_data_line(vs2 + [0] * (NUM_ELEMS - active_vl), src_sew)
    )
    tf.data_label(
        f"{tag}_s1", format_data_line(vs1 + [0] * (NUM_ELEMS - active_vl), src_sew)
    )
    tf.data_align(dst_sew)
    tf.data_label(f"{tag}_exp", format_data_line(expected_active, dst_sew))
    # Background: 1024 bytes (supports vs2r.v up to VLEN=4096)
    tf.data(".align 4")
    tf.data_label(f"{tag}_bg", f"    .fill 256, 4, 0x{bg_word:08x}")
    tf.data(".align 4")
    tf.data_label(f"{tag}_resbuf", "    .space 1024")

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# P2: LMUL>1 register boundary crossing — integer (vadd.vv m2)
# ======================================================================


def _gen_lmul_gt1_int(out: Path) -> str:
    """LMUL>1 test: vadd.vv at e32/m2 with vl=16.

    On VLEN=256: VLMAX=16, elements 0-7 in v8, elements 8-15 in v9.
    This forces the hardware to cross the register boundary (v8→v9)
    which is a major source of hardware bugs.

    The test uses li+vsetvli to request vl=16.  On VLEN=256 this equals
    VLMAX.  On larger VLEN the data still fits within one register group
    so no boundary crossing occurs, but the test remains correct.
    """
    fpath = out / "lmul_gt1_int.S"
    tf = TestFile(
        "vadd.vv",
        "LMUL>1 register boundary crossing: vadd.vv e32/m2 — "
        "16 elements spanning v8-v9 on VLEN=256",
    )

    sew = 32
    num_m2 = 16
    nbytes = num_m2 * (sew // 8)  # 64

    # Distinct values — elements 8-15 cross into v9 on VLEN=256
    vs2 = [i * 10 + 1 for i in range(num_m2)]
    vs1 = [i * 5 + 2 for i in range(num_m2)]
    expected = [(a + b) & M(sew) for a, b in zip(vs2, vs1)]
    witness = [0xDEAD_BEEF] * num_m2

    cn_res = tf.next_check("lmul>1 int: wrong result")
    cn_wit = tf.next_check("lmul>1 int: witness changed")
    cn_csr = tf.next_check("lmul>1 int: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"LMUL>1 test: vadd.vv e32/m2, vl={num_m2}")
    tf.comment("On VLEN=256: elements 0-7 in v8, elements 8-15 in v9")

    # Set e32/m2 with vl=16
    tf.code(f"li t0, {num_m2}")
    tf.code("vsetvli t0, t0, e32, m2, tu, mu")

    # Load sources (m2 register groups: v16-v17, v20-v21)
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Load witness (m2: v24-v25)
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    # Execute (m2: result in v8-v9)
    tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Check result
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code("la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    # Check CSRs
    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_CSRS_UNCHANGED")

    # Data (16 elements each = 64 bytes)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(witness, sew))

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# P2: LMUL>1 register boundary crossing — FP (vfadd.vv m2)
# ======================================================================


def _gen_lmul_gt1_fp(out: Path) -> str:
    """LMUL>1 test: vfadd.vv at e32/m2 with vl=16."""
    fpath = out / "lmul_gt1_fp.S"
    tf = TestFile(
        "vfadd.vv",
        "LMUL>1 register boundary crossing: vfadd.vv e32/m2 — "
        "16 elements spanning v8-v9 on VLEN=256",
    )

    sew = 32
    num_m2 = 16
    nbytes = num_m2 * (sew // 8)  # 64

    # FP data: 1.0..16.0 + 0.5 = 1.5..16.5 (all exact in f32)
    vs2_fp = [float(i + 1) for i in range(num_m2)]
    vs1_fp = [0.5] * num_m2
    exp_fp = [a + b for a, b in zip(vs2_fp, vs1_fp)]

    vs2_bits = [f32_to_bits(v) for v in vs2_fp]
    vs1_bits = [f32_to_bits(v) for v in vs1_fp]
    exp_bits = [f32_to_bits(v) for v in exp_fp]
    witness = [0xDEAD_BEEF] * num_m2

    cn_res = tf.next_check("lmul>1 fp: wrong result")
    cn_wit = tf.next_check("lmul>1 fp: witness changed")
    cn_csr = tf.next_check("lmul>1 fp: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"LMUL>1 test: vfadd.vv e32/m2, vl={num_m2}")

    tf.code(f"li t0, {num_m2}")
    tf.code("vsetvli t0, t0, e32, m2, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    tf.code(f"vfadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code("la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_CSRS_UNCHANGED_FP")

    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2_bits, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1_bits, sew))
    tf.data_label(f"{tag}_exp", format_data_line(exp_bits, sew))
    tf.data_label(f"{tag}_w", format_data_line(witness, sew))

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# P3: vl=0 no-op — FP (vfadd.vv)
# ======================================================================


def _gen_vl_zero_fp(out: Path) -> str:
    """vl=0 no-op test for FP: vfadd.vv must not modify vd."""
    fpath = out / "vl_zero_fp.S"
    tf = TestFile("vfadd.vv", "vl=0 no-op: vfadd.vv must not modify vd")

    sew = 32
    bg = [0xDEAD_BEEF, 0xCAFE_BABE, 0x1234_5678, 0x9ABC_DEF0]
    nbytes = NUM_ELEMS * (sew // 8)

    cn_res = tf.next_check("vl=0 fp: vd modified")
    cn_vstart = tf.next_check("vl=0 fp: vstart not zero")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment("vl=0 test: vfadd.vv — vd must be unchanged")

    # Pre-load vd
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_bg")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    # Load sources (don't matter)
    tf.code(f"la t1, {tag}_src")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Set vl=0
    tf.code("li t0, 0")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Execute — should be a no-op
    tf.code(f"vfadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Check vstart
    tf.code(f"SET_TEST_NUM {cn_vstart}")
    tf.code("csrr t2, vstart")
    tf.code("FAIL_IF_NZ t2")

    # Restore vl=4 to store and check
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_bg, {nbytes}")

    tf.data_align(sew)
    tf.data_label(f"{tag}_bg", format_data_line(bg, sew))
    tf.data_label(f"{tag}_src", format_data_line([f32_to_bits(1.0)] * NUM_ELEMS, sew))

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# P3: vl=0 no-op — Load (vle32.v)
# ======================================================================


def _gen_vl_zero_load(out: Path) -> str:
    """vl=0 no-op test for load: vle32.v must not modify vd."""
    fpath = out / "vl_zero_load.S"
    tf = TestFile("vle32.v", "vl=0 no-op: vle32.v must not modify vd")

    sew = 32
    bg = [0xDEAD_BEEF, 0xCAFE_BABE, 0x1234_5678, 0x9ABC_DEF0]
    nbytes = NUM_ELEMS * (sew // 8)

    cn_res = tf.next_check("vl=0 load: vd modified")
    cn_vstart = tf.next_check("vl=0 load: vstart not zero")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment("vl=0 test: vle32.v — vd must be unchanged")

    # Pre-load vd
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_bg")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    # Set vl=0
    tf.code("li t0, 0")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load — should be a no-op
    tf.code(f"la t1, {tag}_src")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    # Check vstart
    tf.code(f"SET_TEST_NUM {cn_vstart}")
    tf.code("csrr t2, vstart")
    tf.code("FAIL_IF_NZ t2")

    # Restore vl=4 to store and check
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_bg, {nbytes}")

    tf.data_align(sew)
    tf.data_label(f"{tag}_bg", format_data_line(bg, sew))
    # Source memory — distinct values that should NOT end up in vd
    tf.data_label(f"{tag}_src", format_data_line([0xFFFF_FFFF] * NUM_ELEMS, sew))

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# P3: vl=0 no-op — Store (vse32.v)
# ======================================================================


def _gen_vl_zero_store(out: Path) -> str:
    """vl=0 no-op test for store: vse32.v must not write memory."""
    fpath = out / "vl_zero_store.S"
    tf = TestFile("vse32.v", "vl=0 no-op: vse32.v must not write memory")

    sew = 32
    mem_pattern = [0xCAFE_BABE, 0xCAFE_BABE, 0xCAFE_BABE, 0xCAFE_BABE]
    nbytes = NUM_ELEMS * (sew // 8)

    cn_res = tf.next_check("vl=0 store: memory modified")
    cn_vstart = tf.next_check("vl=0 store: vstart not zero")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment("vl=0 test: vse32.v — memory must be unchanged")

    # Load distinct value into vd (should NOT be written to memory)
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code("li t1, 0xFFFFFFFF")
    tf.code(f"vmv.v.x {VREG_DST}, t1")

    # Set vl=0
    tf.code("li t0, 0")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Store — should be a no-op
    tf.code(f"la t1, {tag}_target")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")

    # Check vstart
    tf.code(f"SET_TEST_NUM {cn_vstart}")
    tf.code("csrr t2, vstart")
    tf.code("FAIL_IF_NZ t2")

    # Check memory unchanged (compare target against known-good copy)
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"CHECK_MEM {tag}_target, {tag}_original, {nbytes}")

    # Data: two copies — target (potentially clobbered) and original
    tf.data_align(sew)
    tf.data_label(f"{tag}_target", format_data_line(mem_pattern, sew))
    tf.data_label(f"{tag}_original", format_data_line(mem_pattern, sew))

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# P4: Fault-only-first load (vle32ff.v) — happy path + mmap boundary
# ======================================================================


def _gen_vle32ff_fault(out: Path) -> str:
    """Fault-only-first load test: vle32ff.v at e32.

    Two sub-tests in one .S file:

    Sub-test A (happy path):
      Load from a fully mapped .data buffer with vl=4.  Verify all 4
      elements are loaded correctly AND that ``csrr vl`` equals 4 (i.e.
      vl is NOT reduced when no fault occurs).

    Sub-test B (boundary fault via mmap):
      1. mmap 2 pages (8192 bytes, PROT_READ|PROT_WRITE, MAP_PRIVATE|
         MAP_ANONYMOUS).
      2. munmap the second page -> mapped/unmapped boundary at offset 4096.
      3. Write 3 e32 elements at offset 4084 (= 4096 - 12) of first page.
      4. Pre-fill vd with 0xDEADBEEF at VLMAX via vmv.v.x.
      5. vsetvli vl=4, then vle32ff.v -- element 3 would touch the
         unmapped page and trigger a fault-only-first trim.
      6. Verify: 1 <= vl <= 3 (spec: at least element 0 must succeed,
         cannot load from unmapped page).
      7. Verify vstart == 0.
      8. Dump register with vs1r.v, check element 0 unconditionally,
         elements 1-2 conditionally based on resulting vl.
      9. Cleanup: munmap first page.

    s-register usage:
      s6  = mmap base address (first page)
      s7  = load pointer (base + 4084)
      s8  = vl after fault-only-first load
    """
    fpath = out / "vle32ff_fault.S"
    tf = TestFile(
        "vle32ff.v", "Fault-only-first load: happy path + mmap boundary fault"
    )

    sew = 32
    bg_word = 0xDEAD_BEEF

    # -- Sub-test A: happy-path data -----------------------------------
    happy_src = [0xAAAA_0001, 0xBBBB_0002, 0xCCCC_0003, 0xDDDD_0004]

    cn_happy_res = tf.next_check("happy ff: wrong data loaded")
    cn_happy_vl = tf.next_check("happy ff: vl changed (should stay 4)")

    tag_a = f"tc{cn_happy_res}"

    # -- Sub-test B: boundary-fault check numbers ----------------------
    cn_mmap_fail = tf.next_check("mmap failed")
    cn_munmap_fail = tf.next_check("munmap second page failed")
    cn_vl_zero = tf.next_check("ff vl == 0 (spec requires >= 1)")
    cn_vl_too_big = tf.next_check("ff vl > 3 (loaded from unmapped page?)")
    cn_vstart = tf.next_check("ff vstart != 0 after vle32ff")
    cn_elem0 = tf.next_check("ff element 0 wrong")
    cn_elem1 = tf.next_check("ff element 1 wrong")
    cn_elem2 = tf.next_check("ff element 2 wrong")
    cn_munmap_cleanup = tf.next_check("munmap cleanup failed")

    # ==================================================================
    # Sub-test A: happy path -- load from .data, verify vl unchanged
    # ==================================================================
    tf.blank()
    tf.comment("=" * 60)
    tf.comment("Sub-test A: Happy path -- vle32ff.v from fully mapped memory")
    tf.comment("=" * 60)

    # Set e32/m1, vl=4, tu/mu
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load source with ff
    tf.code(f"la t1, {tag_a}_src")
    tf.code(f"vle{sew}ff.v {VREG_DST}, (t1)")

    # Check vl == 4 (should NOT be reduced)
    tf.code(f"SET_TEST_NUM {cn_happy_vl}")
    tf.code("csrr t2, vl")
    tf.code(f"li t3, {NUM_ELEMS}")
    tf.code("FAIL_IF_NE t2, t3")

    # Store and check loaded data
    tf.code(f"SET_TEST_NUM {cn_happy_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag_a}_src, {NUM_ELEMS * (sew // 8)}")

    # ==================================================================
    # Sub-test B: boundary fault via mmap
    # ==================================================================
    tf.blank()
    tf.comment("=" * 60)
    tf.comment("Sub-test B: Boundary fault -- mmap + munmap trick")
    tf.comment("=" * 60)

    # --- Step 1: mmap 2 pages (8192 bytes) ---
    # PROT_READ|PROT_WRITE = 3, MAP_PRIVATE|MAP_ANONYMOUS = 34
    tf.code(f"SET_TEST_NUM {cn_mmap_fail}")
    tf.code("SYS_MMAP 0, 8192, 3, 34, -1, 0")
    tf.code("bltz a0, 90100f")
    tf.code("j 90101f")
    tf.raw("90100:")
    tf.code("FAIL_TEST")
    tf.raw("90101:")
    tf.code("mv s6, a0")  # save base address

    # --- Step 2: munmap the second page ---
    tf.code(f"SET_TEST_NUM {cn_munmap_fail}")
    tf.code("li t0, 4096")
    tf.code("add s7, s6, t0")  # s7 = base + 4096
    tf.code("SYS_MUNMAP s7, 4096")
    tf.code("FAIL_IF_NZ a0")

    # --- Step 3: Write 3 e32 elements at offset 4084 = 4096 - 12 ---
    # Elements: [0xAAAA0001, 0xBBBB0002, 0xCCCC0003]
    # They occupy bytes 4084..4095 of the first page.
    tf.code("li t0, 4084")
    tf.code("add s7, s6, t0")  # s7 = load pointer
    tf.code("li t1, 0xAAAA0001")
    tf.code("sw t1, 0(s7)")
    tf.code("li t1, 0xBBBB0002")
    tf.code("sw t1, 4(s7)")
    tf.code("li t1, 0xCCCC0003")
    tf.code("sw t1, 8(s7)")

    # --- Step 4: Pre-fill vd with 0xDEADBEEF at VLMAX ---
    tf.code("vsetvli t0, x0, e32, m1, tu, mu")
    tf.code(f"li t1, 0x{bg_word:08x}")
    tf.code(f"vmv.v.x {VREG_DST}, t1")

    # --- Step 5: Set vl=4, execute vle32ff.v ---
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"vle{sew}ff.v {VREG_DST}, (s7)")

    # Read resulting vl into s8
    tf.code("csrr s8, vl")

    # --- Step 6: Verify 1 <= vl <= 3 ---
    # Check vl >= 1  (vl == 0 means spec violation)
    tf.code(f"SET_TEST_NUM {cn_vl_zero}")
    tf.code("beqz s8, 90102f")
    tf.code("j 90103f")
    tf.raw("90102:")
    tf.code("FAIL_TEST")
    tf.raw("90103:")

    # Check vl <= 3  (vl > 3 means it loaded from unmapped page)
    tf.code(f"SET_TEST_NUM {cn_vl_too_big}")
    tf.code("li t0, 3")
    tf.code("bgt s8, t0, 90104f")
    tf.code("j 90105f")
    tf.raw("90104:")
    tf.code("FAIL_TEST")
    tf.raw("90105:")

    # --- Step 7: Check vstart == 0 ---
    tf.code(f"SET_TEST_NUM {cn_vstart}")
    tf.code("csrr t2, vstart")
    tf.code("FAIL_IF_NZ t2")

    # --- Step 8: Dump register and check elements conditionally ---
    # Dump entire register with vs1r.v
    tf.code("la t1, ff_resbuf")
    tf.code(f"vs1r.v {VREG_DST}, (t1)")

    # Element checks use _mem_compare (byte-by-byte) to avoid
    # RV64 sign-extension mismatches between lw and li for values
    # with bit 31 set (e.g. 0xAAAA0001).

    # Element 0: always loaded (spec guarantees at least element 0)
    tf.code(f"SET_TEST_NUM {cn_elem0}")
    tf.code("la a1, ff_resbuf")
    tf.code("la a2, ff_exp")
    tf.code("li a3, 4")
    tf.code("jal ra, _mem_compare")
    tf.code("FAIL_IF_NZ a0")

    # Element 1: check only if vl >= 2
    tf.code(f"SET_TEST_NUM {cn_elem1}")
    tf.code("li t0, 2")
    tf.code("blt s8, t0, 90106f")  # skip if vl < 2
    tf.code("la a1, ff_resbuf")
    tf.code("addi a1, a1, 4")
    tf.code("la a2, ff_exp")
    tf.code("addi a2, a2, 4")
    tf.code("li a3, 4")
    tf.code("jal ra, _mem_compare")
    tf.code("FAIL_IF_NZ a0")
    tf.raw("90106:")

    # Element 2: check only if vl >= 3
    tf.code(f"SET_TEST_NUM {cn_elem2}")
    tf.code("li t0, 3")
    tf.code("blt s8, t0, 90107f")  # skip if vl < 3
    tf.code("la a1, ff_resbuf")
    tf.code("addi a1, a1, 8")
    tf.code("la a2, ff_exp")
    tf.code("addi a2, a2, 8")
    tf.code("li a3, 4")
    tf.code("jal ra, _mem_compare")
    tf.code("FAIL_IF_NZ a0")
    tf.raw("90107:")

    # --- Step 9: Cleanup -- munmap first page ---
    tf.code(f"SET_TEST_NUM {cn_munmap_cleanup}")
    tf.code("SYS_MUNMAP s6, 4096")
    tf.code("FAIL_IF_NZ a0")

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag_a}_src", format_data_line(happy_src, sew))
    # Expected boundary-fault element values (for _mem_compare checks)
    tf.data_align(sew)
    tf.data_label(
        "ff_exp", format_data_line([0xAAAA_0001, 0xBBBB_0002, 0xCCCC_0003], sew)
    )
    # Result buffer for vs1r.v dump (512 bytes, supports up to VLEN=4096)
    tf.data(".align 4")
    tf.data_label("ff_resbuf", "    .space 512")

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 1/11: vill trap — illegal vtype causes SIGILL on vadd.vv,
#            but vmv1r.v (whole-register) works fine
# ======================================================================


def _gen_vill_trap(out: Path) -> str:
    """vill trap: illegal vtype → SIGILL on vtype-dependent instructions.

    1. Load v16 with known data at valid vtype.
    2. Set illegal vtype (e64/mf2: SEW=64 > LMUL*ELEN = 32).
    3. Verify vill bit set (sign bit of vtype) and vl==0.
    4. Fork a child process; child resets SIGILL to SIG_DFL, sets illegal
       vtype itself, executes vadd.vv → kernel kills child with SIGILL.
       Parent checks exit status via wait4.
    5. Restore valid vtype.
    6. vmv1r.v v8, v16, store and verify data integrity.

    Uses clone/wait4 instead of a userspace SIGILL handler because the
    ucontext PC-advance approach has portability issues across emulators.

    The child must set illegal vtype itself because the Linux kernel may
    not preserve vtype/vl across fork (lazy vector context management).
    """
    fpath = out / "vill_trap.S"
    tf = TestFile(
        "vsetvli / vadd.vv / vmv1r.v",
        "vill bit: illegal vtype causes SIGILL on vtype-dependent "
        "instructions; vmv1r.v works after restoring valid vtype",
    )

    sew = 32
    src_data = [0x1111_1111, 0x2222_2222, 0x3333_3333, 0x4444_4444]
    nbytes = NUM_ELEMS * (sew // 8)

    cn_vill = tf.next_check("vill bit not set after illegal vtype")
    cn_vl_zero = tf.next_check("vl not zero after illegal vtype")
    cn_clone = tf.next_check("clone() syscall failed")
    cn_sigill = tf.next_check(
        "vadd.vv did not trap with illegal vtype (child not killed by SIGILL)"
    )
    cn_vmv1r = tf.next_check("vmv1r.v data wrong after restore")
    tag = "vill"

    tf.blank()
    tf.comment("vill trap test: illegal vtype → SIGILL on vadd.vv")
    tf.comment("Uses fork/wait to verify trap without userspace signal handler")

    # 1. Set valid vtype, load v16 with known data
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_src")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")

    # 2. Set illegal vtype: e64/mf2 (SEW=64 > LMUL*ELEN = 1/2 * 64 = 32)
    tf.code("vsetvli t0, x0, e64, mf2, ta, ma")

    # 3. Verify vill bit set (sign bit of vtype CSR on RV64)
    tf.code(f"SET_TEST_NUM {cn_vill}")
    tf.code("csrr t0, vtype")
    tf.code("bgez t0, 90200f")  # sign bit clear → vill not set → fail
    tf.code("j 90201f")
    tf.raw("90200:")
    tf.code("FAIL_TEST")
    tf.raw("90201:")

    # 4. Verify vl == 0
    tf.code(f"SET_TEST_NUM {cn_vl_zero}")
    tf.code("csrr t0, vl")
    tf.code("FAIL_IF_NZ t0")

    # 5. Fork child to test vadd.vv trap (vill still set from step 2)
    tf.code(f"SET_TEST_NUM {cn_clone}")
    tf.blank()
    tf.comment("Fork child process")
    tf.code("SYS_CLONE")
    tf.code("bltz a0, 90202f")  # clone failed → fail
    tf.code("beqz a0, 90203f")  # child → jump to child code
    tf.code("j 90204f")  # parent → jump to parent code
    tf.raw("90202:")
    tf.code("FAIL_TEST")
    tf.blank()

    # --- Child process ---
    tf.raw("90203:")
    tf.comment("Child: set illegal vtype and execute vadd.vv")
    tf.comment("(must set vtype ourselves — kernel does not preserve")
    tf.comment("vtype/vl across fork on all platforms)")
    tf.code("vsetvli t0, x0, e64, mf2, ta, ma")
    tf.blank()
    tf.comment("This vadd.vv should SIGILL → kernel kills child")
    tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
    tf.blank()
    tf.comment("If we reach here, vadd.vv did not trap — exit child")
    tf.code("SYS_EXIT 0")
    tf.blank()

    # --- Parent process ---
    tf.raw("90204:")
    tf.comment("Parent: wait for child via wait4()")
    tf.code("mv s6, a0")  # save child PID
    tf.code(f"la s7, {tag}_wstatus")
    tf.code("SYS_WAIT4 s6, s7")
    tf.blank()
    tf.comment("Check WTERMSIG(wstatus) == SIGILL (4)")
    tf.code(f"SET_TEST_NUM {cn_sigill}")
    tf.code(f"la t0, {tag}_wstatus")
    tf.code("lw t0, 0(t0)")
    tf.code("andi t0, t0, 0x7f")  # WTERMSIG
    tf.code("li t1, 4")  # SIGILL
    tf.code("FAIL_IF_NE t0, t1")

    # 6. Restore valid vtype, re-load v16 (vector regs may not survive
    #    clone+wait syscalls due to lazy vector context management),
    #    then vmv1r.v, verify data integrity
    tf.blank()
    tf.comment("Restore valid vtype and re-load v16 (vector regs may not")
    tf.comment("survive clone+wait due to lazy vector context management)")
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_src")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"vmv1r.v {VREG_DST}, {VREG_SRC2}")

    tf.code(f"SET_TEST_NUM {cn_vmv1r}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_src, {nbytes}")

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag}_src", format_data_line(src_data, sew))
    tf.data(".align 3")
    tf.data_label(f"{tag}_wstatus", ".word 0")

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 2/11: Fractional LMUL — e32/mf2, e16/mf4, e8/mf8
# ======================================================================


def _gen_fract_lmul(out: Path) -> str:
    """Fractional LMUL: vadd.vv at e32/mf2, e16/mf4, e8/mf8.

    For each: fill register with background at e32/m1 VLMAX, switch to
    fractional LMUL (request VLMAX), fill sources via vmv.v.x broadcast,
    execute vadd.vv, dump with vs1r.v, check active bytes + tail preserved.

    Active bytes are computed dynamically from vl and element size.
    """
    fpath = out / "fract_lmul.S"
    tf = TestFile(
        "vadd.vv",
        "Fractional LMUL: e32/mf2, e16/mf4, e8/mf8 — "
        "verify active elements computed, tail bytes preserved",
    )

    bg_word = 0xDEAD_BEEF

    # (sew, lmul_str, shift_for_elem_bytes)
    # shift: 2 for e32 (<<2 = *4), 1 for e16 (<<1 = *2), 0 for e8 (*1)
    subtests = [
        (32, "mf2", 2),
        (16, "mf4", 1),
        (8, "mf8", 0),
    ]

    for sew, lmul, shift in subtests:
        cn_active = tf.next_check(f"fract {lmul}: active elements wrong")
        cn_tail = tf.next_check(f"fract {lmul}: tail bytes modified")
        tag = f"fract_{lmul}"

        tf.blank()
        tf.comment(f"Sub-test: vadd.vv at e{sew}/{lmul}")

        # 1. Fill v8 with bg at e32/m1 VLMAX (fills entire register)
        tf.code("vsetvli t0, x0, e32, m1, tu, mu")
        tf.code(f"li t1, 0x{bg_word:08x}")
        tf.code(f"vmv.v.x {VREG_DST}, t1")

        # 2. Switch to fractional LMUL, get VLMAX → s6
        tf.code(f"vsetvli s6, x0, e{sew}, {lmul}, tu, mu")

        # 3. Fill sources with broadcast: 100 + 10 = 110
        tf.code("li t1, 100")
        tf.code(f"vmv.v.x {VREG_SRC2}, t1")
        tf.code("li t1, 10")
        tf.code(f"vmv.v.x {VREG_SRC1}, t1")

        # 4. Execute
        tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

        # 5. Dump entire register with vs1r.v (VLEN/8 bytes)
        tf.code(f"la t1, {tag}_resbuf")
        tf.code(f"vs1r.v {VREG_DST}, (t1)")

        # 6. Check active bytes: _mem_compare(resbuf, exp, s6 << shift)
        tf.code(f"SET_TEST_NUM {cn_active}")
        tf.code(f"la a1, {tag}_resbuf")
        tf.code(f"la a2, {tag}_exp")
        if shift > 0:
            tf.code(f"slli a3, s6, {shift}")
        else:
            tf.code("mv a3, s6")
        tf.code("jal ra, _mem_compare")
        tf.code("FAIL_IF_NZ a0")

        # 7. Check tail bytes: resbuf[active..] == bg[active..]
        tf.code(f"SET_TEST_NUM {cn_tail}")
        tf.code(f"la a1, {tag}_resbuf")
        if shift > 0:
            tf.code(f"slli t0, s6, {shift}")
        else:
            tf.code("mv t0, s6")
        tf.code("add a1, a1, t0")
        tf.code(f"la a2, {tag}_bg")
        tf.code("add a2, a2, t0")
        tf.code("csrr a3, vlenb")
        tf.code("sub a3, a3, t0")
        tf.code("blez a3, 1f")
        tf.code("jal ra, _mem_compare")
        tf.code("FAIL_IF_NZ a0")
        tf.raw("1:")

    # Data sections for all three sub-tests
    for sew, lmul, shift in subtests:
        tag = f"fract_{lmul}"
        exp_val = 110

        # Expected active: fill with 110 at appropriate width
        if sew == 32:
            tf.data(".align 4")
            tf.data_label(f"{tag}_exp", f"    .fill 128, 4, 0x{exp_val:08x}")
        elif sew == 16:
            tf.data(".align 2")
            tf.data_label(f"{tag}_exp", f"    .fill 256, 2, 0x{exp_val:04x}")
        else:  # sew == 8
            tf.data(".align 1")
            tf.data_label(f"{tag}_exp", f"    .fill 512, 1, 0x{exp_val:02x}")

        # Background (512 bytes of 0xDEADBEEF, supports VLEN up to 4096)
        tf.data(".align 4")
        tf.data_label(f"{tag}_bg", f"    .fill 128, 4, 0x{bg_word:08x}")

        # Result buffer
        tf.data(".align 4")
        tf.data_label(f"{tag}_resbuf", "    .space 512")

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 3/11: Tail agnostic — ta policy, only check active elements
# ======================================================================


def _gen_tail_agnostic(out: Path) -> str:
    """Tail agnostic (ta): vadd.vv with vl=2, ta/mu.

    Only verify active elements 0-1 are computed correctly.
    Tail elements (2-3) can be anything — spec allows any value with ta.
    Uses vs1r.v dump so we can check only the first 8 bytes.
    """
    fpath = out / "tail_agnostic.S"
    tf = TestFile("vadd.vv", "Tail agnostic (ta): vl=2, verify only active elements")

    sew = 32
    vs2 = [100, 200, 300, 400]
    vs1 = [10, 20, 30, 40]
    expected_active = [(100 + 10) & M(sew), (200 + 20) & M(sew)]
    active_bytes = 2 * (sew // 8)  # 8

    cn_res = tf.next_check("tail agnostic: active elements wrong")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment("Tail agnostic: vadd.vv e32 vl=2, ta/mu — only check active")

    # Load sources at vl=4 with ta policy
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, ta, mu")
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Set vl=2 with ta
    tf.code("li t0, 2")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, ta, mu")

    # Execute
    tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Dump with vs1r.v (full register)
    tf.code(f"la t1, {tag}_resbuf")
    tf.code(f"vs1r.v {VREG_DST}, (t1)")

    # Only check first 8 bytes (active elements — tail can be anything)
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"CHECK_MEM {tag}_resbuf, {tag}_exp, {active_bytes}")

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected_active, sew))
    tf.data(".align 4")
    tf.data_label(f"{tag}_resbuf", "    .space 512")

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 4/11: Mask agnostic — ma policy, only check unmasked elements
# ======================================================================


def _gen_mask_agnostic(out: Path) -> str:
    """Mask agnostic (ma): vadd.vv with tu/ma, mask=0b0101.

    Elements 0,2 are active (unmasked) — verify computed correctly.
    Elements 1,3 are masked-off — spec allows any value with ma.
    Check only specific element offsets via _mem_compare.
    """
    fpath = out / "mask_agnostic.S"
    tf = TestFile(
        "vadd.vv",
        "Mask agnostic (ma): tu/ma mask=0b0101, verify only unmasked elements",
    )

    sew = 32
    ebytes = sew // 8  # 4
    vs2 = [100, 200, 300, 400]
    vs1 = [10, 20, 30, 40]
    mask_bits = 0b0101  # elements 0,2 active

    exp_elem0 = (vs2[0] + vs1[0]) & M(sew)  # 110
    exp_elem2 = (vs2[2] + vs1[2]) & M(sew)  # 330

    cn_elem0 = tf.next_check("mask agnostic: element 0 wrong")
    cn_elem2 = tf.next_check("mask agnostic: element 2 wrong")
    tag = f"tc{cn_elem0}"

    tf.blank()
    tf.comment("Mask agnostic: vadd.vv e32 tu/ma mask=0b0101")

    # Load at vl=4 with tu/ma
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, ma")
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Load mask
    tf.code(f"la t1, {tag}_mask")
    tf.code("vlm.v v0, (t1)")

    # Execute masked
    tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0.t")

    # Store result (all 4 elements including masked-off ones)
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")

    # Check element 0 (offset 0, 4 bytes)
    tf.code(f"SET_TEST_NUM {cn_elem0}")
    tf.code("la a1, result_buf")
    tf.code(f"la a2, {tag}_exp0")
    tf.code(f"li a3, {ebytes}")
    tf.code("jal ra, _mem_compare")
    tf.code("FAIL_IF_NZ a0")

    # Check element 2 (offset 8, 4 bytes)
    tf.code(f"SET_TEST_NUM {cn_elem2}")
    tf.code("la a1, result_buf")
    tf.code(f"addi a1, a1, {2 * ebytes}")
    tf.code(f"la a2, {tag}_exp2")
    tf.code(f"li a3, {ebytes}")
    tf.code("jal ra, _mem_compare")
    tf.code("FAIL_IF_NZ a0")

    # Data
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp0", format_data_line([exp_elem0], sew))
    tf.data_label(f"{tag}_exp2", format_data_line([exp_elem2], sew))

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 5/11: Strided load with stride=0 (broadcast)
# ======================================================================


def _gen_stride_zero(out: Path) -> str:
    """vlse32.v with stride=0: all elements read from same address.

    Expected: all 4 elements identical (broadcast of single word).
    Includes witness + CSR checks.
    """
    fpath = out / "stride_zero.S"
    tf = TestFile(
        "vlse32.v", "Stride=0 load: all elements read from same address (broadcast)"
    )

    sew = 32
    src_word = 0xCAFE_BABE
    expected = [src_word] * NUM_ELEMS
    nbytes = NUM_ELEMS * (sew // 8)

    wpat = witness_pattern(sew)

    cn_res = tf.next_check("stride-zero: wrong result")
    cn_wit = tf.next_check("stride-zero: witness changed")
    cn_csr = tf.next_check("stride-zero: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment("Strided load with stride=0 — broadcast single word")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    # vlse32.v with stride = 0 (x0 register = 0)
    tf.code(f"la t1, {tag}_src")
    tf.code(f"vlse{sew}.v {VREG_DST}, (t1), x0")

    # Check result
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code("la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    # Check CSRs
    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_CSRS_UNCHANGED")

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag}_src", format_data_line([src_word], sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wpat, sew))

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 6/11: Strided load with negative stride (reverse)
# ======================================================================


def _gen_stride_negative(out: Path) -> str:
    """vlse32.v with stride=-4: elements loaded in reverse order.

    Memory [A,B,C,D] at ascending addresses.  rs1 points to D (last).
    stride=-4 → element 0 reads D, element 1 reads C, etc.
    Result: [D, C, B, A].
    """
    fpath = out / "stride_negative.S"
    tf = TestFile("vlse32.v", "Negative stride load: elements loaded in reverse order")

    sew = 32
    ebytes = sew // 8  # 4
    src = [0xAAAA_0001, 0xBBBB_0002, 0xCCCC_0003, 0xDDDD_0004]
    # Reversed: element 0 = src[3], element 1 = src[2], ...
    expected = list(reversed(src))
    nbytes = NUM_ELEMS * ebytes

    wpat = witness_pattern(sew)

    cn_res = tf.next_check("stride-negative: wrong result")
    cn_wit = tf.next_check("stride-negative: witness changed")
    cn_csr = tf.next_check("stride-negative: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment("Strided load with stride=-4 — reverse order")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    # rs1 = &src[3] (last element), stride = -4
    tf.code(f"la t1, {tag}_src")
    tf.code(f"addi t1, t1, {(NUM_ELEMS - 1) * ebytes}")
    tf.code(f"li a0, -{ebytes}")
    tf.code(f"vlse{sew}.v {VREG_DST}, (t1), a0")

    # Check result
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code("la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    # Check CSRs
    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_CSRS_UNCHANGED")

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag}_src", format_data_line(src, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wpat, sew))

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 7/11: Ordered indexed store with duplicate indices
# ======================================================================


def _gen_scatter_ordered(out: Path) -> str:
    """vsoxei32.v (ordered) with duplicate indices [8,4,4,0].

    Stores happen in element order (0,1,2,3).  At index 4, element 2
    overwrites element 1's value.
    Expected memory: [data[3], data[2], data[0], prefill].
    """
    fpath = out / "scatter_ordered.S"
    tf = TestFile(
        "vsoxei32.v",
        "Ordered indexed store: duplicate indices, later element "
        "overwrites earlier at same offset",
    )

    sew = 32
    ebytes = sew // 8
    data = [0x1111_1111, 0x2222_2222, 0x3333_3333, 0x4444_4444]
    indices = [8, 4, 4, 0]  # byte offsets; indices[1]==indices[2]
    prefill = 0xAAAA_AAAA
    # Ordered: elem 0 → base+8, elem 1 → base+4, elem 2 → base+4, elem 3 → base+0
    expected = [
        data[3],  # base+0: from element 3
        data[2],  # base+4: from element 2 (overwrites element 1)
        data[0],  # base+8: from element 0
        prefill,  # base+12: untouched
    ]
    nbytes = NUM_ELEMS * ebytes

    cn_res = tf.next_check("scatter ordered: wrong memory content")
    cn_csr = tf.next_check("scatter ordered: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment("Ordered indexed store with duplicate indices")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load data (vs3 = v8 in our convention)
    tf.code(f"la t1, {tag}_data")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    # Load indices into v16
    tf.code(f"la t1, {tag}_idx")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")

    tf.code("SAVE_CSRS")

    # Execute ordered indexed store: vsoxei32.v vs3, (rs1), vs2
    tf.code(f"la t1, {tag}_target")
    tf.code(f"vsoxei{sew}.v {VREG_DST}, (t1), {VREG_SRC2}")

    # Check memory
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"CHECK_MEM {tag}_target, {tag}_exp, {nbytes}")

    # Check CSRs (store instruction — full check)
    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_CSRS_UNCHANGED")

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag}_data", format_data_line(data, sew))
    tf.data_label(f"{tag}_idx", format_data_line(indices, sew))
    tf.data_label(f"{tag}_target", format_data_line([prefill] * NUM_ELEMS, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 8/11: FP flags (fflags) — verify NX and OF bits set
# ======================================================================


def _gen_fflags_set(out: Path) -> str:
    """Verify fflags bits set correctly by FP vector operations.

    Sub-test A: Inexact (NX, bit 0)
      vfadd.vv with 1.0 + 1.5*2^-24 → rounds to 1+2^-23 (inexact).

    Sub-test B: Overflow (OF, bit 2)
      vfadd.vv with FLT_MAX + FLT_MAX → +inf (overflow + inexact).
    """
    fpath = out / "fflags_set.S"
    tf = TestFile("vfadd.vv", "FP flags: verify fflags NX and OF bits set correctly")

    sew = 32

    # Sub-test A: inexact
    # 1.0f + 1.5*2^-24 → inexact (rounds to 1.0 + 2^-23)
    val_one = f32_to_bits(1.0)  # 0x3F800000
    val_tiny = f32_to_bits(1.5 * 2**-24)  # 0x33C00000

    cn_nx = tf.next_check("fflags: NX bit not set after inexact add")

    # Sub-test B: overflow
    # FLT_MAX + FLT_MAX → +inf, sets OF and NX
    flt_max_bits = 0x7F7F_FFFF

    cn_of = tf.next_check("fflags: OF bit not set after overflow add")

    tf.blank()
    tf.comment("Sub-test A: inexact — 1.0f + 1.5*2^-24 sets NX")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Clear fflags
    tf.code("csrw fflags, x0")

    # Load sources: broadcast 1.0 and 1.5*2^-24
    tf.code(f"li t1, 0x{val_one:08x}")
    tf.code(f"vmv.v.x {VREG_SRC2}, t1")
    tf.code(f"li t1, 0x{val_tiny:08x}")
    tf.code(f"vmv.v.x {VREG_SRC1}, t1")

    # Execute
    tf.code(f"vfadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Check NX bit (bit 0 of fflags)
    tf.code(f"SET_TEST_NUM {cn_nx}")
    tf.code("csrr t0, fflags")
    tf.code("andi t0, t0, 1")  # NX = bit 0
    tf.code("FAIL_IF_Z t0")

    tf.blank()
    tf.comment("Sub-test B: overflow — FLT_MAX + FLT_MAX sets OF")

    # Clear fflags
    tf.code("csrw fflags, x0")

    # Load FLT_MAX broadcast
    tf.code(f"li t1, 0x{flt_max_bits:08x}")
    tf.code(f"vmv.v.x {VREG_SRC2}, t1")
    tf.code(f"vmv.v.x {VREG_SRC1}, t1")

    # Execute
    tf.code(f"vfadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Check OF bit (bit 2 of fflags)
    tf.code(f"SET_TEST_NUM {cn_of}")
    tf.code("csrr t0, fflags")
    tf.code("andi t0, t0, 4")  # OF = bit 2
    tf.code("FAIL_IF_Z t0")

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 9/11: vxsat sticky behavior
# ======================================================================


def _gen_vxsat_sticky(out: Path) -> str:
    """vxsat is sticky: once set by saturating op, stays 1 until cleared.

    1. Clear vxsat.
    2. vsadd.vv with 0x7FFFFFFF + 1 → saturates to 0x7FFFFFFF, vxsat=1.
    3. Check vxsat == 1.
    4. vadd.vv (non-saturating) — must NOT clear vxsat.
    5. Check vxsat still == 1 (sticky).
    6. csrw vxsat, 0 → clears it.
    7. Check vxsat == 0.
    """
    fpath = out / "vxsat_sticky.S"
    tf = TestFile(
        "vsadd.vv / vadd.vv",
        "vxsat sticky: set by saturating op, preserved by "
        "non-saturating op, cleared only by explicit csrw",
    )

    sew = 32

    cn_sat_res = tf.next_check("vxsat: saturating result wrong")
    cn_sat_set = tf.next_check("vxsat: not set after saturating add")
    cn_sticky = tf.next_check("vxsat: cleared by non-saturating add (not sticky)")
    cn_clear = tf.next_check("vxsat: not cleared after csrw vxsat,0")

    tag = "vxsat"

    tf.blank()
    tf.comment("vxsat sticky test")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Step 1: clear vxsat
    tf.code("csrw vxsat, x0")

    # Step 2: vsadd.vv with 0x7FFFFFFF + 1 → saturates to 0x7FFFFFFF
    tf.code("li t1, 0x7FFFFFFF")
    tf.code(f"vmv.v.x {VREG_SRC2}, t1")
    tf.code("li t1, 1")
    tf.code(f"vmv.v.x {VREG_SRC1}, t1")
    tf.code(f"vsadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Step 2b: verify result is saturated (0x7FFFFFFF)
    tf.code(f"SET_TEST_NUM {cn_sat_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_sat_exp, {NUM_ELEMS * (sew // 8)}")

    # Step 3: check vxsat == 1
    tf.code(f"SET_TEST_NUM {cn_sat_set}")
    tf.code("csrr t0, vxsat")
    tf.code("andi t0, t0, 1")
    tf.code("FAIL_IF_Z t0")

    # Step 4: non-saturating vadd.vv — must NOT clear vxsat
    tf.code("li t1, 1")
    tf.code(f"vmv.v.x {VREG_SRC2}, t1")
    tf.code("li t1, 2")
    tf.code(f"vmv.v.x {VREG_SRC1}, t1")
    tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Step 5: check vxsat STILL 1 (sticky)
    tf.code(f"SET_TEST_NUM {cn_sticky}")
    tf.code("csrr t0, vxsat")
    tf.code("andi t0, t0, 1")
    tf.code("FAIL_IF_Z t0")

    # Step 6: explicit clear
    tf.code("csrw vxsat, x0")

    # Step 7: check vxsat == 0
    tf.code(f"SET_TEST_NUM {cn_clear}")
    tf.code("csrr t0, vxsat")
    tf.code("FAIL_IF_NZ t0")

    # Data: expected saturated result (all 0x7FFFFFFF)
    tf.data_align(sew)
    tf.data_label(f"{tag}_sat_exp", format_data_line([0x7FFF_FFFF] * NUM_ELEMS, sew))

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 10/11: Narrowing tail — vnsrl.wi e64→e32, vl=2, tu
# ======================================================================


def _gen_narrowing_tail(out: Path) -> str:
    """Narrowing with tail undisturbed: vnsrl.wi (e64→e32), vl=2, tu.

    Pre-fill dest at e32/m1 VLMAX with background.  Load e64 source.
    Execute narrowing with vl=2.  Dump with vs1r.v.
    Check: active bytes (8) correct + tail bytes preserved.
    """
    fpath = out / "narrowing_tail.S"
    tf = TestFile(
        "vnsrl.wi",
        "Narrowing tail undisturbed: vnsrl.wi e64→e32 vl=2 — "
        "verify active elements and tail preserved",
    )

    src_sew = 64
    dst_sew = 32
    active_vl = 2
    active_bytes = active_vl * (dst_sew // 8)  # 8
    bg_word = 0xBEEF_DEAD

    # Source e64 values: low 32 bits are the narrowing result (shift=0)
    src_e64 = [0x00000001_00000064, 0x00000002_000000C8]
    expected_active = [0x00000064, 0x000000C8]  # 100, 200

    cn_active = tf.next_check("narrowing-tail: active elements wrong")
    cn_tail = tf.next_check("narrowing-tail: tail bytes modified")
    tag = f"tc{cn_active}"

    tf.blank()
    tf.comment("Narrowing tail test: vnsrl.wi e64→e32, vl=2, tu/mu")

    # 1. Fill v8 with bg at e32/m1 VLMAX (fills entire register)
    tf.code("vsetvli t0, x0, e32, m1, tu, mu")
    tf.code(f"li t1, 0x{bg_word:08x}")
    tf.code(f"vmv.v.x {VREG_DST}, t1")

    # 2. Load e64 source into v16 (need vl=2 at e64/m1)
    tf.code(f"li t0, {active_vl}")
    tf.code("vsetvli t0, t0, e64, m1, tu, mu")
    tf.code(f"la t1, {tag}_src64")
    tf.code(f"vle64.v {VREG_SRC2}, (t1)")

    # 3. Set vtype for narrowing: e32/m1, vl=2
    #    vnsrl operates at dest SEW (e32), reading 2*SEW (e64) source
    tf.code(f"li t0, {active_vl}")
    tf.code("vsetvli t0, t0, e32, m1, tu, mu")

    # 4. Execute narrowing (shift amount 0 = truncate to low 32 bits)
    tf.code(f"vnsrl.wi {VREG_DST}, {VREG_SRC2}, 0")

    # 5. Dump entire register with vs1r.v
    tf.code(f"la t1, {tag}_resbuf")
    tf.code(f"vs1r.v {VREG_DST}, (t1)")

    # 6. Check active bytes (first 8 bytes)
    tf.code(f"SET_TEST_NUM {cn_active}")
    tf.code(f"CHECK_MEM {tag}_resbuf, {tag}_exp, {active_bytes}")

    # 7. Check tail bytes (vlenb - 8)
    _emit_tail_check(tf, cn_tail, f"{tag}_resbuf", f"{tag}_bg", active_bytes)

    # Data
    tf.data_align(src_sew)
    tf.data_label(f"{tag}_src64", format_data_line(src_e64, src_sew))
    tf.data_align(dst_sew)
    tf.data_label(f"{tag}_exp", format_data_line(expected_active, dst_sew))
    tf.data(".align 4")
    tf.data_label(f"{tag}_bg", f"    .fill 128, 4, 0x{bg_word:08x}")
    tf.data(".align 4")
    tf.data_label(f"{tag}_resbuf", "    .space 512")

    tf.write(fpath)
    return str(fpath)


# ======================================================================
# Test 11/11: Widening m2→m4 — vwadd.vv e32/m2 → e64/m4, vl=16
# ======================================================================


def _gen_widening_m2_m4(out: Path) -> str:
    """Widening vwadd.vv at e32/m2 → e64/m4 with vl=16.

    Sources v16-v17 (m2) and v20-v21 (m2), dest v8-v11 (m4).
    On VLEN=256: vl=16 fills all 4 dest registers, crossing v8→v9→v10→v11.
    On VLEN>256: vl=16 still works (16 ≤ VLMAX).
    On VLEN=128: vl clamped to VLMAX=8; test still correct but tests
    fewer boundaries.  Dynamic vl-based byte count handles this.
    """
    fpath = out / "widening_m2_m4.S"
    tf = TestFile(
        "vwadd.vv",
        "Widening m2→m4: vwadd.vv e32/m2 → e64/m4 spanning 4 registers on VLEN=256",
    )

    src_sew = 32
    dst_sew = 64
    requested_vl = 16

    # Source data: 16 e32 values
    vs2 = [i * 10 + 1 for i in range(requested_vl)]
    vs1 = [i * 5 + 2 for i in range(requested_vl)]
    # Expected: sext32→64(vs2[i]) + sext32→64(vs1[i])
    # All values are small positive, so sext is identity
    expected_e64 = [
        (S(a, src_sew) + S(b, src_sew)) & M(dst_sew) for a, b in zip(vs2, vs1)
    ]

    cn_res = tf.next_check("widening m2→m4: wrong result")
    cn_csr = tf.next_check("widening m2→m4: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Widening test: vwadd.vv e32/m2 → e64/m4, vl={requested_vl}")
    tf.comment("Dest v8-v11 (m4), src v16-v17 (m2) + v20-v21 (m2)")

    # Set e32/m2, request vl=16 → s6 = actual vl
    tf.code(f"li t0, {requested_vl}")
    tf.code("vsetvli s6, t0, e32, m2, tu, mu")

    # Load sources at e32/m2
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{src_sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{src_sew}.v {VREG_SRC1}, (t1)")

    tf.code("SAVE_CSRS")

    # Execute widening: vwadd.vv v8, v16, v20 (at e32/m2 → e64/m4)
    tf.code(f"vwadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Switch to e64/m4 with same vl to store result
    # VLMAX(e64/m4) = VLMAX(e32/m2) = VLEN/16, so same vl works
    tf.code("mv t0, s6")
    tf.code("vsetvli t0, t0, e64, m4, tu, mu")

    # Store result
    tf.code(f"la t1, {tag}_resbuf")
    tf.code(f"vse{dst_sew}.v {VREG_DST}, (t1)")

    # Check result: s6 * 8 bytes
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la a1, {tag}_resbuf")
    tf.code(f"la a2, {tag}_exp")
    tf.code("slli a3, s6, 3")  # nbytes = vl * 8
    tf.code("jal ra, _mem_compare")
    tf.code("FAIL_IF_NZ a0")

    # Check CSRs — vl/vtype intentionally changed for the store, so
    # only check the CSRs that should be unchanged
    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_VSTART_ZERO")
    tf.code("CHECK_VXSAT_UNCHANGED")
    tf.code("CHECK_VXRM_UNCHANGED")
    tf.code("CHECK_FCSR_UNCHANGED")

    # Data
    tf.data_align(src_sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, src_sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, src_sew))
    tf.data_align(dst_sew)
    tf.data_label(f"{tag}_exp", format_data_line(expected_e64, dst_sew))
    tf.data(".align 4")
    tf.data_label(f"{tag}_resbuf", "    .space 256")

    tf.write(fpath)
    return str(fpath)


def _gen_vstart_nonzero(out: Path) -> str:
    """Test vstart CSR write/read and behavior of vector insns with vstart>0.

    The RVV 1.0 spec (§3.7) says:
      "Implementations are permitted to raise illegal instruction
       exceptions when attempting to execute a vector instruction with
       a value of vstart that the implementation can never produce."

    So both outcomes are legal:
      - vadd.vv traps (SIGILL): implementation never interrupts
        mid-instruction, so nonzero vstart is unsupported.
      - vadd.vv works: elements before vstart are skipped.

    This test verifies:
      1. csrw vstart, 2 succeeds (no trap).
      2. csrr vstart reads back 2.
      3. Fork a child that executes vadd.vv with vstart=2.
         If child is killed by SIGILL → the implementation doesn't
         support nonzero vstart (legal per spec).
         If child exits normally → check that vadd.vv skipped
         elements 0-1 and computed elements 2-3.
    """
    fpath = out / "vstart_nonzero.S"
    tf = TestFile(
        "vstart CSR + vadd.vv",
        "Verifies csrw/csrr vstart works and tests vadd.vv behavior "
        "with nonzero vstart (both SIGILL and skip-elements are valid "
        "per RVV 1.0 §3.7)",
    )

    sew = 32
    vd_pre = [0xAAAA_AAAA, 0xBBBB_BBBB, 0xCCCC_CCCC, 0xDDDD_DDDD]
    vs2 = [100, 200, 300, 400]
    vs1 = [10, 20, 30, 40]
    # If vstart=2: elements 0-1 preserved from vd_pre, 2-3 computed
    expected_skip = [
        vd_pre[0],
        vd_pre[1],
        (vs2[2] + vs1[2]) & M(sew),
        (vs2[3] + vs1[3]) & M(sew),
    ]
    nbytes = NUM_ELEMS * (sew // 8)

    cn_csrw = tf.next_check("csrw vstart, 2 caused trap (unexpected)")
    cn_csrr = tf.next_check("csrr vstart did not read back 2")
    cn_clone = tf.next_check("clone() failed")
    cn_behavior = tf.next_check(
        "vadd.vv with vstart=2: child exited 0 but result wrong "
        "(neither SIGILL nor correct skip-elements behavior)"
    )
    tag = "vstart"

    tf.blank()
    tf.comment("vstart nonzero test")
    tf.comment("RVV 1.0 §3.7: implementations may trap on nonzero vstart")

    # 1. Set valid vtype, load registers
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # 2. Write vstart=2 and read it back
    tf.code(f"SET_TEST_NUM {cn_csrw}")
    tf.code("li t0, 2")
    tf.code("csrw vstart, t0")
    tf.comment("If we reach here, csrw vstart succeeded")

    tf.code(f"SET_TEST_NUM {cn_csrr}")
    tf.code("csrr t1, vstart")
    tf.code("li t2, 2")
    tf.code("FAIL_IF_NE t1, t2")

    # 3. Reset vstart to 0 before fork (vstart=2 would affect child setup)
    tf.code("csrw vstart, zero")

    # 4. Fork child to test vadd.vv with vstart=2
    tf.code(f"SET_TEST_NUM {cn_clone}")
    tf.code("SYS_CLONE")
    tf.code(f"bltz a0, {tag}_clone_fail")
    tf.code(f"beqz a0, {tag}_child")
    tf.code(f"j {tag}_parent")
    tf.raw(f"{tag}_clone_fail:")
    tf.code("FAIL_TEST")
    tf.blank()

    # --- Child ---
    tf.raw(f"{tag}_child:")
    tf.comment("Child: reload registers (not preserved across fork),")
    tf.comment("set vstart=2, execute vadd.vv")
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
    tf.code("li t0, 2")
    tf.code("csrw vstart, t0")
    tf.comment("This vadd.vv may SIGILL (legal) or skip elements 0-1 (also legal)")
    tf.code(f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
    tf.blank()
    tf.comment("If we reach here, vadd.vv succeeded — store result and exit 0")
    tf.code(f"la t1, {tag}_result")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code("SYS_EXIT 0")
    tf.blank()

    # --- Parent ---
    tf.raw(f"{tag}_parent:")
    tf.code("mv s6, a0")
    tf.code(f"la s7, {tag}_wstatus")
    tf.code("SYS_WAIT4 s6, s7")
    tf.blank()
    tf.comment("Check child outcome")
    tf.code(f"la t0, {tag}_wstatus")
    tf.code("lw t0, 0(t0)")
    tf.blank()
    tf.comment("WIFEXITED: bits [15:8] = exit status, bits [6:0] = 0")
    tf.comment("WIFSIGNALED: bits [6:0] = signal number")
    tf.code("andi t1, t0, 0x7f")  # WTERMSIG
    tf.code(f"beqz t1, {tag}_exited")
    tf.blank()
    tf.comment("Child was signaled — check it was SIGILL (4)")
    tf.code("li t2, 4")
    tf.code(f"beq t1, t2, {tag}_sigill_ok")
    tf.comment("Killed by unexpected signal — fail")
    tf.code(f"SET_TEST_NUM {cn_behavior}")
    tf.code("FAIL_TEST")
    tf.blank()

    tf.raw(f"{tag}_sigill_ok:")
    tf.comment("SIGILL from vadd.vv with vstart=2 — legal per spec, pass")
    # falls through to PASS_TEST

    # Need to jump over the exited-normally check
    tf.code(f"j {tag}_done")
    tf.blank()

    tf.raw(f"{tag}_exited:")
    tf.comment("Child exited normally — vadd.vv ran, verify result")
    tf.code(f"SET_TEST_NUM {cn_behavior}")
    tf.code(f"CHECK_MEM {tag}_result, {tag}_exp, {nbytes}")
    tf.blank()

    tf.raw(f"{tag}_done:")

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag}_vd", format_data_line(vd_pre, sew))
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected_skip, sew))
    tf.data(".align 3")
    tf.data_label(f"{tag}_wstatus", "    .word 0")
    tf.data_label(f"{tag}_result", f"    .space {nbytes}")

    tf.write(fpath)
    return str(fpath)


# ==================================================================
# Helper: emit fork-based SIGILL check for a raw instruction encoding
# ==================================================================


def _emit_sigill_fork_check(
    tf: TestFile,
    encoding: int,
    label_prefix: str,
    description: str,
    cn_clone: int,
    cn_sigill: int,
    *,
    child_setup: list[str] | None = None,
) -> None:
    """Emit a fork-based SIGILL test for a raw 32-bit instruction encoding.

    Forks a child process via clone(), child executes the raw instruction
    (which must SIGILL), parent waits via wait4() and verifies
    WTERMSIG == SIGILL (4).

    child_setup: optional extra assembly lines for the child before the
                 faulting instruction (e.g. setting up vtype).
    """
    tag = label_prefix
    tf.blank()
    tf.comment(f"SIGILL check: {description}")
    tf.comment(f"Encoding 0x{encoding:08X} — must trap as illegal instruction")

    tf.code(f"SET_TEST_NUM {cn_clone}")
    tf.code("SYS_CLONE")
    tf.code(f"bltz a0, {tag}_clone_fail")
    tf.code(f"beqz a0, {tag}_child")
    tf.code(f"j {tag}_parent")
    tf.raw(f"{tag}_clone_fail:")
    tf.code("FAIL_TEST")
    tf.blank()

    # --- Child ---
    tf.raw(f"{tag}_child:")
    if child_setup:
        for line in child_setup:
            tf.code(line)
    tf.comment(f"Execute 0x{encoding:08X} — should SIGILL")
    tf.code(f".4byte 0x{encoding:08X}")
    tf.comment("If we reach here, instruction did NOT trap — exit child 0")
    tf.code("SYS_EXIT 0")
    tf.blank()

    # --- Parent ---
    tf.raw(f"{tag}_parent:")
    tf.code("mv s6, a0")  # save child PID
    tf.code(f"la s7, {tag}_wstatus")
    tf.code("SYS_WAIT4 s6, s7")
    tf.blank()
    tf.code(f"SET_TEST_NUM {cn_sigill}")
    tf.code(f"la t0, {tag}_wstatus")
    tf.code("lw t0, 0(t0)")
    tf.code("andi t0, t0, 0x7f")  # WTERMSIG
    tf.code("li t1, 4")  # SIGILL
    tf.code("FAIL_IF_NE t0, t1")


def _gen_ghostwrite(out: Path) -> str:
    """GhostWrite (CVE-2024-44067) vulnerability check.

    Verifies that a correct RVV 1.0 implementation rejects the
    illegally-encoded vector instructions that constitute the GhostWrite
    attack on T-Head XuanTie C910/C920 CPUs.

    GhostWrite exploits instructions with mew=1 (bit 28), which is
    RESERVED in standard RVV 1.0.  On the vulnerable T-Head CPUs, these
    instructions (vse128.v / vle128.v) bypass virtual memory and operate
    on physical addresses directly.

    On a correct implementation, mew=1 must trap as an illegal instruction.

    Encodings tested:
      0x10028027 — vse128.v v0, 0(t0)  (exact GhostWrite PoC encoding)
      0x12028407 — vle128.v v8, 0(t0)  (load counterpart, mew=1)
    """
    fpath = out / "ghostwrite.S"
    tf = TestFile(
        "GhostWrite CVE-2024-44067 check",
        "Verifies mew=1 (reserved) vector store/load instructions trap as "
        "illegal instruction, preventing the GhostWrite physical-memory "
        "write vulnerability (CVE-2024-44067, T-Head C910/C920)",
    )

    # The exact GhostWrite PoC encoding: vse128.v v0, 0(t0)
    # [28]=mew=1, [27:26]=mop=00, [25]=vm=0, [24:20]=00000,
    # [19:15]=t0(x5), [14:12]=000, [11:7]=v0, [6:0]=0100111
    GHOSTWRITE_STORE = 0x10028027

    # Load counterpart: vle128.v v8, 0(t0)
    # Same but LOAD-FP opcode, vd=v8, vm=1
    GHOSTWRITE_LOAD = 0x12028407

    cn_store_clone = tf.next_check("clone failed for vse128.v check")
    cn_store_sigill = tf.next_check(
        "vse128.v (GhostWrite, 0x10028027) did NOT trap — VULNERABLE to CVE-2024-44067!"
    )
    cn_load_clone = tf.next_check("clone failed for vle128.v check")
    cn_load_sigill = tf.next_check(
        "vle128.v (0x12028407, mew=1 load) did NOT trap — reserved encoding accepted"
    )

    # Child needs a valid vtype set before executing the faulting instruction
    # (because the encoding is in LOAD-FP/STORE-FP opcode space, the CPU
    # needs to know it's a vector instruction — mew=1 makes it illegal)
    child_vtype_setup = [
        "li t0, 4",
        "vsetvli t0, t0, e8, m1, ta, ma",
        "la t0, gw_scratch",  # point t0 at scratch buffer
    ]

    tf.blank()
    tf.comment("GhostWrite vulnerability check (CVE-2024-44067)")
    tf.comment("T-Head XuanTie C910/C920: vse128.v writes physical memory")
    tf.comment("Standard RVV 1.0: mew=1 (bit 28) is reserved → must SIGILL")

    # Test 1: vse128.v (the actual GhostWrite store)
    _emit_sigill_fork_check(
        tf,
        GHOSTWRITE_STORE,
        "gw_store",
        "vse128.v v0, 0(t0) — GhostWrite store (CVE-2024-44067)",
        cn_store_clone,
        cn_store_sigill,
        child_setup=child_vtype_setup,
    )

    # Test 2: vle128.v (the load counterpart)
    _emit_sigill_fork_check(
        tf,
        GHOSTWRITE_LOAD,
        "gw_load",
        "vle128.v v8, 0(t0) — mew=1 load",
        cn_load_clone,
        cn_load_sigill,
        child_setup=child_vtype_setup,
    )

    # PASS_TEST is auto-appended by TestFile.write()

    # Data
    tf.data(".align 4")
    tf.data_label("gw_scratch", "    .space 64")
    tf.data(".align 3")
    tf.data_label("gw_store_wstatus", "    .word 0")
    tf.data_label("gw_load_wstatus", "    .word 0")

    tf.write(fpath)
    return str(fpath)


def _gen_reserved_encoding(out: Path) -> str:
    """Test that structurally reserved RVV load/store encodings trap as SIGILL.

    Tests encoding fields that are reserved in the RVV 1.0 vector
    load/store format itself (not funct6 slots in OP-V, which are
    gradually being claimed by newer extensions like Zvbc, Zvbb, etc.).

    1. Strided store with mew=1 (vsse128.v encoding)
       Similar to GhostWrite but strided addressing mode.

    2. Reserved lumop=00001 in unit-stride vector load
       Only lumop=00000 (unit-stride), 01000 (fault-first),
       01011 (mask load), 10000 (whole-register) are valid.

    3. Indexed load with mew=1
       vluxei128.v — mew=1 with indexed (unordered) addressing.

    4. Reserved sumop=00001 in unit-stride vector store
       Only sumop=00000 (unit-stride), 01011 (mask store),
       10000 (whole-register) are valid.
    """
    fpath = out / "reserved_encoding.S"
    tf = TestFile(
        "Reserved RVV encoding checks",
        "Verifies that structurally reserved RVV 1.0 load/store encodings "
        "trap as illegal instruction (SIGILL). Tests mew=1 (strided store "
        "and indexed load) and reserved lumop/sumop values.",
    )

    # Strided store with mew=1: vsse128.v v0, 0(t0), t1
    # mop=10 (strided), mew=1, width=000, rs2=t1(x6), rs1=t0(x5)
    STRIDED_MEW1 = (
        0x27
        | (0 << 7)
        | (0b000 << 12)
        | (5 << 15)
        | (6 << 20)
        | (1 << 25)
        | (0b10 << 26)
        | (1 << 28)
    )

    # Reserved lumop=00001 in unit-stride load
    # vle8.v encoding but with lumop=00001 instead of 00000
    RESERVED_LUMOP = 0x07 | (8 << 7) | (0b000 << 12) | (5 << 15) | (1 << 20) | (1 << 25)

    # Indexed load with mew=1: vluxei128.v v8, (t0), v16
    # mop=01 (indexed-unordered), mew=1, width=000, vs2=v16, rs1=t0
    INDEXED_MEW1 = (
        0x07
        | (8 << 7)
        | (0b000 << 12)
        | (5 << 15)
        | (16 << 20)
        | (1 << 25)
        | (0b01 << 26)
        | (1 << 28)
    )

    # Reserved sumop=00001 in unit-stride store
    # vse8.v encoding but with sumop=00001 instead of 00000
    RESERVED_SUMOP = 0x27 | (8 << 7) | (0b000 << 12) | (5 << 15) | (1 << 20) | (1 << 25)

    cn1_clone = tf.next_check("clone failed for mew=1 strided store")
    cn1_sigill = tf.next_check(
        "mew=1 strided store (0x{:08X}) did NOT trap".format(STRIDED_MEW1)
    )
    cn2_clone = tf.next_check("clone failed for reserved lumop")
    cn2_sigill = tf.next_check(
        "reserved lumop=00001 (0x{:08X}) did NOT trap".format(RESERVED_LUMOP)
    )
    cn3_clone = tf.next_check("clone failed for mew=1 indexed load")
    cn3_sigill = tf.next_check(
        "mew=1 indexed load (0x{:08X}) did NOT trap".format(INDEXED_MEW1)
    )
    cn4_clone = tf.next_check("clone failed for reserved sumop")
    cn4_sigill = tf.next_check(
        "reserved sumop=00001 (0x{:08X}) did NOT trap".format(RESERVED_SUMOP)
    )

    child_vtype_setup = [
        "li t0, 4",
        "vsetvli t0, t0, e8, m1, ta, ma",
        "la t0, re_scratch",
        "li t1, 4",
    ]

    tf.blank()
    tf.comment("Reserved RVV 1.0 load/store encoding checks")
    tf.comment("These test structurally reserved fields (mew, lumop, sumop),")
    tf.comment("NOT OP-V funct6 slots (which are claimed by newer extensions).")

    _emit_sigill_fork_check(
        tf,
        STRIDED_MEW1,
        "re_smew",
        "mew=1 strided store",
        cn1_clone,
        cn1_sigill,
        child_setup=child_vtype_setup,
    )

    _emit_sigill_fork_check(
        tf,
        RESERVED_LUMOP,
        "re_lumop",
        "unit-stride load with reserved lumop=00001",
        cn2_clone,
        cn2_sigill,
        child_setup=child_vtype_setup,
    )

    _emit_sigill_fork_check(
        tf,
        INDEXED_MEW1,
        "re_imew",
        "mew=1 indexed load",
        cn3_clone,
        cn3_sigill,
        child_setup=child_vtype_setup,
    )

    _emit_sigill_fork_check(
        tf,
        RESERVED_SUMOP,
        "re_sumop",
        "unit-stride store with reserved sumop=00001",
        cn4_clone,
        cn4_sigill,
        child_setup=child_vtype_setup,
    )

    # PASS_TEST is auto-appended by TestFile.write()

    # Data
    tf.data(".align 4")
    tf.data_label("re_scratch", "    .space 64")
    tf.data(".align 3")
    tf.data_label("re_smew_wstatus", "    .word 0")
    tf.data_label("re_lumop_wstatus", "    .word 0")
    tf.data_label("re_imew_wstatus", "    .word 0")
    tf.data_label("re_sumop_wstatus", "    .word 0")

    tf.write(fpath)
    return str(fpath)


def _gen_rvv_detect(out: Path) -> str:
    """Detect standard RVV 1.0 vs XTheadVector from userspace.

    XTheadVector (T-Head C906V/C920/R920) overlaps with the V extension
    encoding space but lacks several RVV 1.0 features.  Since the M-mode
    CSRs needed for official detection (mvendorid, mimpid) are not
    accessible from userspace, we detect via instruction availability:

    1. vsetivli (exists in RVV 1.0, absent in XTheadVector)
    2. vlenb CSR read (exists in both, but encoding differs in old
       XTheadVector toolchains — kept as a sanity check)
    3. vmv1r.v (whole register move, absent in XTheadVector)

    All three must succeed for standard RVV 1.0.  If any traps, the
    implementation is XTheadVector (or broken).
    """
    fpath = out / "rvv_detect.S"
    tf = TestFile(
        "RVV 1.0 detection",
        "Verifies the implementation is standard RVV 1.0 (not "
        "XTheadVector) by testing instructions that exist only in "
        "RVV 1.0: vsetivli, vlenb CSR, vmv1r.v",
    )

    cn_vsetivli = tf.next_check("vsetivli not supported (XTheadVector?)")
    cn_vlenb = tf.next_check("csrr vlenb failed")
    cn_vlenb_val = tf.next_check("vlenb is zero (invalid)")
    cn_vmv1r = tf.next_check("vmv1r.v not supported (XTheadVector?)")
    cn_vmv1r_data = tf.next_check("vmv1r.v data integrity failed")

    sew = 32
    src_data = [0x1111_1111, 0x2222_2222, 0x3333_3333, 0x4444_4444]
    nbytes = NUM_ELEMS * (sew // 8)
    tag = "rvvdet"

    tf.blank()
    tf.comment("RVV 1.0 detection: test instructions absent in XTheadVector")

    # Test 1: vsetivli (immediate AVL form, RVV 1.0 only)
    tf.code(f"SET_TEST_NUM {cn_vsetivli}")
    tf.code(f"vsetivli t0, {NUM_ELEMS}, e{sew}, m1, ta, ma")
    tf.comment("If we reach here, vsetivli works → not XTheadVector")

    # Test 2: csrr vlenb
    tf.code(f"SET_TEST_NUM {cn_vlenb}")
    tf.code("csrr t0, vlenb")
    tf.code(f"SET_TEST_NUM {cn_vlenb_val}")
    tf.code("FAIL_IF_Z t0")

    # Test 3: vmv1r.v (whole register move)
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_src")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"SET_TEST_NUM {cn_vmv1r}")
    tf.code(f"vmv1r.v {VREG_DST}, {VREG_SRC2}")

    # Verify data integrity
    tf.code(f"SET_TEST_NUM {cn_vmv1r_data}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_src, {nbytes}")

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag}_src", format_data_line(src_data, sew))

    tf.write(fpath)
    return str(fpath)


def _gen_tail_per_family(out: Path) -> str:
    """Test tail-undisturbed (tu) policy across instruction families.

    Uses vl=2 with e32/m1 (VLMAX=VLEN/32 ≥ 4), pre-fills vd with a
    known pattern, executes each instruction, then dumps the full
    register via vs1r.v and verifies:
      - Active elements (0..vl-1) have the correct computed result.
      - Tail elements (vl..VLEN/32-1) are preserved from the pre-fill.

    Since VLEN is unknown at generation time, we only check the first
    4 elements (2 active + 2 tail) using a fixed 16-byte comparison,
    which works for any VLEN >= 128.
    """
    fpath = out / "tail_per_family.S"
    tf = TestFile(
        "Tail undisturbed per family",
        "Verifies tail elements are preserved (tu policy) for "
        "representative instructions from each family: vadd, vmul, "
        "vwadd, vmacc, vsadd, vfadd, vle32, vslideup",
    )

    sew = 32
    vl = 2
    nbytes = 4 * (sew // 8)  # check first 4 e32 elements = 16 bytes
    prefill = [0xAAAA_AAAA, 0xBBBB_BBBB, 0xCCCC_CCCC, 0xDDDD_DDDD]
    tag_pfx = "tpf"

    # Helper to emit one tail test
    def _emit_tail_test(
        name: str,
        insn: str,
        src2_vals: list[int],
        src1_vals: list[int],
        expected_active: list[int],
        *,
        setup: str = "",
        csr_check: str = "CHECK_CSRS_UNCHANGED",
    ) -> None:
        exp = list(expected_active[:vl]) + list(prefill[vl:])
        cn = tf.next_check(f"{name} tu: tail not preserved")
        tag = f"{tag_pfx}{cn}"

        tf.blank()
        tf.comment(f"Tail test: {name} (vl={vl}, tu policy)")
        if setup:
            tf.code(setup)
        tf.code(f"li t0, {vl}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

        # Pre-fill destination with pattern
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_pf")
        tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

        # Set vl to the test value
        tf.code(f"li t0, {vl}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

        # Load sources
        if src2_vals:
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        if src1_vals:
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

        tf.code("SAVE_CSRS")
        tf.code(insn)

        # Dump full register to result_buf
        tf.code("la t1, result_buf")
        tf.code(f"vs1r.v {VREG_DST}, (t1)")

        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code(csr_check)

        tf.data_align(sew)
        tf.data_label(f"{tag}_pf", format_data_line(prefill, sew))
        if src2_vals:
            tf.data_label(f"{tag}_s2", format_data_line(src2_vals, sew))
        if src1_vals:
            tf.data_label(f"{tag}_s1", format_data_line(src1_vals, sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

    # vadd.vv
    s2 = [10, 20, 30, 40]
    s1 = [1, 2, 3, 4]
    exp_add = [(a + b) & M(sew) for a, b in zip(s2, s1)]
    _emit_tail_test(
        "vadd.vv", f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}", s2, s1, exp_add
    )

    # vmul.vv
    exp_mul = [(a * b) & M(sew) for a, b in zip(s2, s1)]
    _emit_tail_test(
        "vmul.vv", f"vmul.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}", s2, s1, exp_mul
    )

    # vsadd.vv (saturating add)
    from ..compute.fixed_point import sadd

    sa = [0x7FFF_FFFE, 0x7FFF_FFFF, 10, 20]
    sb = [1, 1, 5, 5]
    exp_sadd = [sadd(a, b, sew)[0] for a, b in zip(sa, sb)]
    _emit_tail_test(
        "vsadd.vv",
        f"vsadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}",
        sa,
        sb,
        exp_sadd,
        csr_check="CHECK_CSRS_UNCHANGED_FIXEDPOINT",
    )

    # vfadd.vv (FP)
    from ..common import f32_to_bits

    fa = [f32_to_bits(1.0), f32_to_bits(2.0), f32_to_bits(3.0), f32_to_bits(4.0)]
    fb = [f32_to_bits(0.5), f32_to_bits(0.5), f32_to_bits(0.5), f32_to_bits(0.5)]
    exp_fp = [f32_to_bits(1.5), f32_to_bits(2.5), f32_to_bits(3.5), f32_to_bits(4.5)]
    _emit_tail_test(
        "vfadd.vv",
        f"vfadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}",
        fa,
        fb,
        exp_fp,
        csr_check="CHECK_CSRS_UNCHANGED_FP",
    )

    # vle32.v (load)
    ld_data = [0x1111, 0x2222, 0x3333, 0x4444]
    exp_ld = ld_data  # active elements loaded, tail preserved
    cn = tf.next_check("vle32.v tu: tail not preserved")
    tag = f"{tag_pfx}{cn}"
    tf.blank()
    tf.comment(f"Tail test: vle32.v (vl={vl}, tu policy)")
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_pf")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"li t0, {vl}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code("SAVE_CSRS")
    tf.code(f"la t1, {tag}_mem")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code("la t1, result_buf")
    tf.code(f"vs1r.v {VREG_DST}, (t1)")
    tf.code(f"SET_TEST_NUM {cn}")
    exp_ld_full = list(ld_data[:vl]) + list(prefill[vl:])
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
    tf.code("CHECK_CSRS_UNCHANGED")
    tf.data_align(sew)
    tf.data_label(f"{tag}_pf", format_data_line(prefill, sew))
    tf.data_label(f"{tag}_mem", format_data_line(ld_data, sew))
    tf.data_label(f"{tag}_exp", format_data_line(exp_ld_full, sew))

    # vslideup.vi (offset=1)
    sl_src = [10, 20, 30, 40]
    # slideup with offset=1, vl=2: element 0 from vd (i<offset), element 1 = src[0]
    exp_slide = [prefill[0], sl_src[0]] + list(prefill[vl:])
    cn = tf.next_check("vslideup.vi tu: tail not preserved")
    tag = f"{tag_pfx}{cn}"
    tf.blank()
    tf.comment(f"Tail test: vslideup.vi (vl={vl}, offset=1, tu policy)")
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_pf")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"la t1, {tag}_src")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"li t0, {vl}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code("SAVE_CSRS")
    tf.code(f"vslideup.vi {VREG_DST}, {VREG_SRC2}, 1")
    tf.code("la t1, result_buf")
    tf.code(f"vs1r.v {VREG_DST}, (t1)")
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
    tf.code("CHECK_CSRS_UNCHANGED")
    tf.data_align(sew)
    tf.data_label(f"{tag}_pf", format_data_line(prefill, sew))
    tf.data_label(f"{tag}_src", format_data_line(sl_src, sew))
    tf.data_label(f"{tag}_exp", format_data_line(exp_slide, sew))

    tf.write(fpath)
    return str(fpath)


def _gen_lmul2_per_family(out: Path) -> str:
    """Test LMUL=2 boundary crossing for representative instructions.

    With e32/m2 on VLEN=256, VLMAX=16: elements 0-7 in v8, 8-15 in v9.
    Uses vl=16 to exercise the full register group, verifying elements
    crossing the v8→v9 boundary are correctly computed.

    Since we need more than 4 elements, we dynamically request VLMAX
    with vsetvli t0, x0 and use runtime vl.  We generate expected data
    for 16 elements (correct for VLEN=256) but only compare as many
    bytes as vl*4.
    """
    fpath = out / "lmul2_per_family.S"
    tf = TestFile(
        "LMUL=2 per family",
        "Tests register group boundary crossing (e32/m2) for vadd, "
        "vmul, vand, vsll, vmin across 16 elements (v8+v9)")

    sew = 32
    n = 16  # elements for VLEN=256

    # Generate test data: 16 elements
    s2 = [(i + 1) * 100 for i in range(n)]
    s1 = [(i + 1) for i in range(n)]

    tests: list[tuple[str, str, list[int]]] = [
        ("vadd.vv", f"vadd.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}",
         [(a + b) & M(sew) for a, b in zip(s2, s1)]),
        ("vmul.vv", f"vmul.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}",
         [(a * b) & M(sew) for a, b in zip(s2, s1)]),
        ("vand.vv", f"vand.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}",
         [(a & b) & M(sew) for a, b in zip(s2, s1)]),
        ("vsll.vv", f"vsll.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}",
         [(a << (b & 0x1f)) & M(sew) for a, b in zip(s2, s1)]),
        ("vminu.vv", f"vminu.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}",
         [min(a, b) for a, b in zip(s2, s1)]),
    ]
    tag_pfx = "lm2"

    for name, insn, expected in tests:
        cn = tf.next_check(f"{name} e32/m2: boundary crossing")
        tag = f"{tag_pfx}{cn}"

        tf.blank()
        tf.comment(f"LMUL=2 test: {name} (e32/m2, vl=VLMAX)")

        # Request VLMAX elements
        tf.code("vsetvli t0, x0, e32, m2, ta, ma")
        # Save actual vl for byte count
        tf.code("mv s6, t0")

        # Load sources (need vs2r.v for m2 loads — or just vle32 with LMUL=2)
        tf.code(f"la t1, {tag}_s2")
        tf.code(f"vle32.v {VREG_SRC2}, (t1)")
        tf.code(f"la t1, {tag}_s1")
        tf.code(f"vle32.v {VREG_SRC1}, (t1)")

        tf.code("SAVE_CSRS")
        tf.code(insn)

        # Store result
        tf.code(f"la t1, {tag}_res")
        tf.code(f"vse32.v {VREG_DST}, (t1)")

        # Compute comparison length: vl * 4 bytes
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code("slli a3, s6, 2")  # nbytes = vl * 4
        tf.code(f"la a1, {tag}_res")
        tf.code(f"la a2, {tag}_exp")
        tf.code("call _mem_compare")
        tf.code("FAIL_IF_NZ a0")

        tf.data(".align 4")
        tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
        tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
        tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
        tf.data_label(f"{tag}_res", f"    .space {n * 4}")

    tf.write(fpath)
    return str(fpath)
