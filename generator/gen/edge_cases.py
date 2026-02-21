from __future__ import annotations

"""Generate tests for architectural edge cases:
- vl=0 no-op behavior (integer, FP, load, store)
- vstart>0 partial execution (commented out - SIGILL on target)
- Tail undisturbed (vl < VLMAX) — basic 4-element check
- Full-VLMAX tail undisturbed — verifies ALL tail bytes via vs1r.v + vlenb
- LMUL>1 register boundary crossing (v8→v9)

These test the vector unit's control logic rather than specific ALU
operations.  We use vadd.vv / vfadd.vv / vle32.v / vse32.v / vwadd.vv
as representative instructions since the behavior is
instruction-independent.
"""

from pathlib import Path

from ..common import NUM_ELEMS, format_data_line, format_hex, M, f32_to_bits
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
    # Test 2: vstart>0 — SKIPPED
    # ==================================================================
    # Writing to vstart from userspace causes SIGILL on this target
    # (Linux kernel traps csrw vstart).  This is a known platform
    # limitation — vstart testing requires bare-metal or kernel support.

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
        vd_pre[1],                      # element 1: inactive (mu), preserved
        vd_pre[2],                      # element 2: tail (tu), preserved
        vd_pre[3],                      # element 3: tail (tu), preserved
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
    tf = TestFile("vadd.vv",
                  "Full-VLMAX tail undisturbed: vadd.vv e32 — "
                  "verifies ALL tail bytes up to VLEN/8")

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
    tf.comment("Fill register with bg at VLMAX, compute 2 elements, "
               "verify tail preserved")

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
    _emit_tail_check(tf, cn_tail, f"{tag}_resbuf", f"{tag}_bg",
                     active_bytes)

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2",
                  format_data_line(vs2 + [0] * (NUM_ELEMS - active_vl), sew))
    tf.data_label(f"{tag}_s1",
                  format_data_line(vs1 + [0] * (NUM_ELEMS - active_vl), sew))
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
    tf = TestFile("vfadd.vv",
                  "Full-VLMAX tail undisturbed: vfadd.vv e32 — "
                  "verifies ALL tail bytes up to VLEN/8")

    sew = 32
    active_vl = 2
    active_bytes = active_vl * (sew // 8)  # 8
    bg_word = 0xDEAD_BEEF

    vs2_fp = [1.0, 2.0]
    vs1_fp = [0.5, 0.25]
    exp_fp = [a + b for a, b in zip(vs2_fp, vs1_fp)]   # [1.5, 2.25]

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
    _emit_tail_check(tf, cn_tail, f"{tag}_resbuf", f"{tag}_bg",
                     active_bytes)

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2",
                  format_data_line(vs2_bits + [0] * (NUM_ELEMS - active_vl),
                                   sew))
    tf.data_label(f"{tag}_s1",
                  format_data_line(vs1_bits + [0] * (NUM_ELEMS - active_vl),
                                   sew))
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
    tf = TestFile("vle32.v",
                  "Full-VLMAX tail undisturbed: vle32.v — "
                  "verifies ALL tail bytes up to VLEN/8")

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
    _emit_tail_check(tf, cn_tail, f"{tag}_resbuf", f"{tag}_bg",
                     active_bytes)

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
    tf = TestFile("vwadd.vv",
                  "Full-VLMAX tail undisturbed: vwadd.vv e32->e64 — "
                  "verifies ALL tail bytes in dest group")

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
    tf.comment(f"Full-VLMAX tail test: vwadd.vv e32→e64, vl={active_vl}, "
               "tu/mu")
    tf.comment("Dest is v8-v9 (EMUL=2). Uses vs2r.v to dump, "
               "2*vlenb for tail.")

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
    _emit_tail_check(tf, cn_tail, f"{tag}_resbuf", f"{tag}_bg",
                     active_bytes, nreg=2)

    # Data
    tf.data_align(src_sew)
    tf.data_label(f"{tag}_s2",
                  format_data_line(vs2 + [0] * (NUM_ELEMS - active_vl),
                                   src_sew))
    tf.data_label(f"{tag}_s1",
                  format_data_line(vs1 + [0] * (NUM_ELEMS - active_vl),
                                   src_sew))
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
    tf = TestFile("vadd.vv",
                  "LMUL>1 register boundary crossing: vadd.vv e32/m2 — "
                  "16 elements spanning v8-v9 on VLEN=256")

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
    tf = TestFile("vfadd.vv",
                  "LMUL>1 register boundary crossing: vfadd.vv e32/m2 — "
                  "16 elements spanning v8-v9 on VLEN=256")

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
    tf.data_label(f"{tag}_src",
                  format_data_line([f32_to_bits(1.0)] * NUM_ELEMS, sew))

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
    tf.data_label(f"{tag}_src",
                  format_data_line([0xFFFF_FFFF] * NUM_ELEMS, sew))

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
    tf = TestFile("vle32ff.v",
                  "Fault-only-first load: happy path + mmap boundary fault")

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
    # mmap(addr=0, length=8192, prot=PROT_READ|PROT_WRITE,
    #       flags=MAP_PRIVATE|MAP_ANONYMOUS, fd=-1, offset=0)
    # Linux RV64: syscall 222
    tf.code(f"SET_TEST_NUM {cn_mmap_fail}")
    tf.code("li a0, 0")                  # addr = NULL (let kernel choose)
    tf.code("li a1, 8192")               # length = 2 pages
    tf.code("li a2, 3")                  # PROT_READ(1) | PROT_WRITE(2)
    tf.code("li a3, 34")                 # MAP_PRIVATE(2) | MAP_ANONYMOUS(32)
    tf.code("li a4, -1")                 # fd = -1
    tf.code("li a5, 0")                  # offset = 0
    tf.code("li a7, 222")               # __NR_mmap
    tf.code("ecall")
    # Check for error (negative a0 on RV64 Linux = error)
    tf.code("bltz a0, 90100f")
    tf.code("j 90101f")
    tf.raw("90100:")
    tf.code("FAIL_TEST")
    tf.raw("90101:")
    # Save base address in s6
    tf.code("mv s6, a0")

    # --- Step 2: munmap the second page ---
    tf.code(f"SET_TEST_NUM {cn_munmap_fail}")
    tf.code("li t0, 4096")
    tf.code("add a0, s6, t0")            # addr = base + 4096
    tf.code("li a1, 4096")               # length = 1 page
    tf.code("li a7, 215")               # __NR_munmap
    tf.code("ecall")
    tf.code("FAIL_IF_NZ a0")

    # --- Step 3: Write 3 e32 elements at offset 4084 = 4096 - 12 ---
    # Elements: [0xAAAA0001, 0xBBBB0002, 0xCCCC0003]
    # They occupy bytes 4084..4095 of the first page.
    tf.code("li t0, 4084")
    tf.code("add s7, s6, t0")            # s7 = load pointer
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
    tf.code("blt s8, t0, 90106f")        # skip if vl < 2
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
    tf.code("blt s8, t0, 90107f")        # skip if vl < 3
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
    tf.code("mv a0, s6")                 # addr = base
    tf.code("li a1, 4096")               # length = 1 page
    tf.code("li a7, 215")               # __NR_munmap
    tf.code("ecall")
    tf.code("FAIL_IF_NZ a0")

    # Data
    tf.data_align(sew)
    tf.data_label(f"{tag_a}_src", format_data_line(happy_src, sew))
    # Expected boundary-fault element values (for _mem_compare checks)
    tf.data_align(sew)
    tf.data_label("ff_exp",
                  format_data_line([0xAAAA_0001, 0xBBBB_0002, 0xCCCC_0003], sew))
    # Result buffer for vs1r.v dump (512 bytes, supports up to VLEN=4096)
    tf.data(".align 4")
    tf.data_label("ff_resbuf", "    .space 512")

    tf.write(fpath)
    return str(fpath)
