from __future__ import annotations

"""Assembly emission helpers – shared code that adds test-case boilerplate
to a :class:`TestFile`.

Every ``emit_*`` function appends code + data for ONE test case (one
combination of SEW + test-vector) and returns the allocated check numbers
so callers can log them.
"""

from .common import (
    NUM_ELEMS, format_data_line, format_hex, witness_pattern,
)
from .testfile import TestFile

# Register allocation (kept constant across all templates):
#   v8       – primary destination
#   v16      – vs2 (first vector source)
#   v20      – vs1 (second vector source)
#   v24      – witness register
#   v0       – mask register (when needed)
#   a0       – scalar operand for VX / VF forms
#   fa0      – FP scalar operand for VF forms

VREG_DST = "v8"
VREG_SRC2 = "v16"
VREG_SRC1 = "v20"
VREG_WITNESS = "v24"


# ===================================================================
# Generic binary: VV form
# ===================================================================

def emit_binop_vv(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    vs1: list[int],
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name}: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name}")

    # Setup vtype
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load sources
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    # Save CSRs
    tf.code("SAVE_CSRS")

    # Execute
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Check result
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    # Check CSRs
    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    # Data
    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Generic binary: VX form
# ===================================================================

def emit_binop_vx(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    scalar: int,
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name}: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name}")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")

    # Load scalar into a0
    tf.code(f"li a0, {format_hex(scalar, sew)}")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, a0")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Generic binary: VI form
# ===================================================================

def emit_binop_vi(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    imm: int,
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name}: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name}")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {imm}")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Widening binary: VV form  (dest is 2*SEW)
# ===================================================================

def emit_widening_vv(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    vs1: list[int],
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    dsew = 2 * sew
    nbytes = NUM_ELEMS * (dsew // 8)
    src_nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew}\u2192e{dsew} {test_name}: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew}\u2192e{dsew} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew}\u2192e{dsew} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name}")

    # Set up source SEW for loading operands
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Load witness (at source SEW)
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    # Execute (instruction uses source SEW, destination EMUL=2)
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Check witness (at source SEW, before vsetvli change)
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {src_nbytes}")

    # Check CSRs (before vsetvli changes vl/vtype)
    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    # Store result at destination width
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    wp = witness_pattern(sew)
    tf.data_align(dsew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, dsew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Compare: VV form  (dest is mask register)
# ===================================================================

def emit_compare_vv(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    vs1: list[int],
    expected_mask: int,
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name}: wrong mask result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name}")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    # Write result to v8 (not v0, to avoid mask-register conflicts)
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Extract mask: store mask byte and check
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vsm.v {VREG_DST}, (t1)")
    tf.code(f"lbu t2, 0(t1)")
    # Only check the low NUM_ELEMS bits
    tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
    tf.code(f"li t3, {expected_mask}")
    tf.code("FAIL_IF_NE t2, t3")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Ternary accumulate: VV form  (vd = f(vd, vs1, vs2))
# ===================================================================

def emit_ternary_vv(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vd_init: list[int],
    vs1: list[int],
    vs2: list[int],
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name}: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name}")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load vd initial value
    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    # Load vs1
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
    # Load vs2
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    # Note: for vmacc/vnmsac, asm is: vmacc.vv vd, vs1, vs2
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC1}, {VREG_SRC2}")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_vd", format_data_line(vd_init, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# FP binary: VV form  (same structure as int binop but with FP CSR check)
# ===================================================================

def emit_fp_binop_vv(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    vs1: list[int],
    expected: list[int],
) -> None:
    emit_binop_vv(
        tf, mnemonic, sew, test_name, vs2, vs1, expected,
        csr_check="CHECK_CSRS_UNCHANGED_FP",
    )


# ===================================================================
# FP binary: VF form
# ===================================================================

def emit_fp_binop_vf(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    scalar_bits: int,
    expected: list[int],
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name}: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name}")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")

    # Load FP scalar
    load_insn = "flw" if sew == 32 else "fld"
    tf.code(f"la t1, {tag}_fsc")
    tf.code(f"{load_insn} fa0, 0(t1)")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, fa0")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_CSRS_UNCHANGED_FP")

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_fsc", format_data_line([scalar_bits], sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Reduction: vs  (scalar result in element 0 of vd)
# ===================================================================

def emit_reduction(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    scalar_init: int,
    vec: list[int],
    expected_scalar: int,
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name}: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name}")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load vector source into vs2
    tf.code(f"la t1, {tag}_vec")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    # Load scalar init into vs1 element 0
    tf.code(f"la t1, {tag}_init")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")
    # Reduction: vd = op(vs1[0], vs2[*]), result in vd[0]
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    # Check element 0 of result
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    # Only check first element
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {sew // 8}")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    # Data
    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_vec", format_data_line(vec, sew))
    # For scalar init, put value in element 0, rest don't matter
    tf.data_label(f"{tag}_init", format_data_line([scalar_init] + [0] * (NUM_ELEMS - 1), sew))
    tf.data_label(f"{tag}_exp", format_data_line([expected_scalar] + [0] * (NUM_ELEMS - 1), sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


def emit_reduction_masked(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    scalar_init: int,
    vec: list[int],
    mask_bits: int,
    expected_scalar: int,
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    """Test masked reduction: only active elements participate."""
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (masked)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_vec")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_init")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
    tf.code(f"la t1, {tag}_mask")
    tf.code("vlm.v v0, (t1)")
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0.t")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code("la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {sew // 8}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code("la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
    tf.data_align(sew)
    tf.data_label(f"{tag}_vec", format_data_line(vec, sew))
    tf.data_label(f"{tag}_init", format_data_line([scalar_init] + [0] * (NUM_ELEMS - 1), sew))
    tf.data_label(f"{tag}_exp", format_data_line([expected_scalar] + [0] * (NUM_ELEMS - 1), sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Mask-register logical  (mm form)
# ===================================================================

def emit_mask_logical(
    tf: TestFile,
    mnemonic: str,
    test_name: str,
    vs2_mask: int,
    vs1_mask: int,
    expected_mask: int,
    vl: int = NUM_ELEMS,
) -> None:
    cn_res = tf.next_check(f"{mnemonic} {test_name}: wrong mask result")
    cn_wit = tf.next_check(f"{mnemonic} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"
    sew = 8  # mask ops use e8
    nbytes = NUM_ELEMS * 1  # e8 → 1 byte per element

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} {test_name}")

    # Use e8/m1 for mask operations
    tf.code(f"li t0, {vl}")
    tf.code(f"vsetvli t0, t0, e8, m1, tu, mu")

    # Load mask sources as byte data
    tf.code(f"la t1, {tag}_m2")
    tf.code(f"vlm.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_m1")
    tf.code(f"vlm.v {VREG_SRC1}, (t1)")

    # Load witness (at e8)
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle8.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vsm.v {VREG_DST}, (t1)")
    tf.code(f"lbu t2, 0(t1)")
    tf.code(f"andi t2, t2, {(1 << vl) - 1}")
    tf.code(f"li t3, {expected_mask & ((1 << vl) - 1)}")
    tf.code("FAIL_IF_NE t2, t3")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse8.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_CSRS_UNCHANGED")

    wp = witness_pattern(sew)
    tf.data(f".align 1")
    tf.data_label(f"{tag}_m2", f"    .byte {vs2_mask}")
    tf.data_label(f"{tag}_m1", f"    .byte {vs1_mask}")
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Compare: VF form  (dest is mask register)
# ===================================================================

def emit_compare_vf(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    scalar_bits: int,
    expected_mask: int,
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name}: wrong mask result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name}")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")

    load_insn = "flw" if sew == 32 else "fld"
    tf.code(f"la t1, {tag}_fsc")
    tf.code(f"{load_insn} fa0, 0(t1)")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, fa0")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vsm.v {VREG_DST}, (t1)")
    tf.code(f"lbu t2, 0(t1)")
    tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
    tf.code(f"li t3, {expected_mask}")
    tf.code("FAIL_IF_NE t2, t3")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_CSRS_UNCHANGED_FP")

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_fsc", format_data_line([scalar_bits], sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Ternary accumulate: VF form  (vd = f(vd, scalar, vs2))
# ===================================================================

def emit_ternary_vf(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vd_init: list[int],
    scalar_bits: int,
    vs2: list[int],
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name}: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name}: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name}: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name}")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load vd initial value
    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    # Load vs2
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")

    # Load FP scalar
    load_insn = "flw" if sew == 32 else "fld"
    tf.code(f"la t1, {tag}_fsc")
    tf.code(f"{load_insn} fa0, 0(t1)")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    # asm: vfmacc.vf vd, rs1, vs2  → vd[i] = rs1 * vs2[i] + vd[i]
    tf.code(f"{mnemonic} {VREG_DST}, fa0, {VREG_SRC2}")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    # Check witness
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_vd", format_data_line(vd_init, sew))
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_fsc", format_data_line([scalar_bits], sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Generic binary: VV overlap form (vd == vs2)
# ===================================================================

def emit_binop_vv_overlap(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    vs1: list[int],
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    """Like emit_binop_vv but uses vd==vs2 (register overlap test)."""
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (vd=vs2)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load vs2 into VREG_DST (overlap: vd == vs2)
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    # Execute with overlap: vd == vs2
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_DST}, {VREG_SRC1}")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Generic binary: VV masked form (v0.t)
# ===================================================================

def emit_binop_vv_masked(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    vs1: list[int],
    vd_init: list[int],
    mask_bits: int,
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    """Test masked operation.  expected[i] = computed if mask bit set, else vd_init[i]."""
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (masked)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # Load sources
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Init vd
    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    # Load mask into v0
    tf.code(f"la t1, {tag}_mask")
    tf.code(f"vlm.v v0, (t1)")

    # Load witness
    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")

    # Execute with mask
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0.t")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_vd", format_data_line(vd_init, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Generic binary: VX overlap form (vd == vs2)
# ===================================================================

def emit_binop_vx_overlap(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    scalar: int,
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    """Like emit_binop_vx but uses vd==vs2 (register overlap test)."""
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (vd=vs2)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"li a0, {format_hex(scalar, sew)}")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_DST}, a0")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Generic binary: VX masked form (v0.t)
# ===================================================================

def emit_binop_vx_masked(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    scalar: int,
    vd_init: list[int],
    mask_bits: int,
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (masked)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"li a0, {format_hex(scalar, sew)}")

    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"la t1, {tag}_mask")
    tf.code(f"vlm.v v0, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, a0, v0.t")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_vd", format_data_line(vd_init, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Generic binary: VI overlap form (vd == vs2)
# ===================================================================

def emit_binop_vi_overlap(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    imm: int,
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (vd=vs2)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_DST}, {imm}")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Generic binary: VI masked form (v0.t)
# ===================================================================

def emit_binop_vi_masked(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    imm: int,
    vd_init: list[int],
    mask_bits: int,
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (masked)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")

    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"la t1, {tag}_mask")
    tf.code(f"vlm.v v0, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {imm}, v0.t")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_vd", format_data_line(vd_init, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# FP binary: VF overlap form (vd == vs2)
# ===================================================================

def emit_fp_binop_vf_overlap(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    scalar_bits: int,
    expected: list[int],
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (vd=vs2)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    load_insn = "flw" if sew == 32 else "fld"
    tf.code(f"la t1, {tag}_fsc")
    tf.code(f"{load_insn} fa0, 0(t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_DST}, fa0")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_CSRS_UNCHANGED_FP")

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_fsc", format_data_line([scalar_bits], sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# FP binary: VF masked form (v0.t)
# ===================================================================

def emit_fp_binop_vf_masked(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    scalar_bits: int,
    vd_init: list[int],
    mask_bits: int,
    expected: list[int],
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (masked)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")

    load_insn = "flw" if sew == 32 else "fld"
    tf.code(f"la t1, {tag}_fsc")
    tf.code(f"{load_insn} fa0, 0(t1)")

    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"la t1, {tag}_mask")
    tf.code(f"vlm.v v0, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, fa0, v0.t")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code("CHECK_CSRS_UNCHANGED_FP")

    wp = witness_pattern(sew)
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_fsc", format_data_line([scalar_bits], sew))
    tf.data_label(f"{tag}_vd", format_data_line(vd_init, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Ternary accumulate: VV overlap form (vd == vs2)
# ===================================================================

def emit_ternary_vv_overlap(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vd_init: list[int],
    vs1: list[int],
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    """Ternary VV with vd==vs2.  asm: vmacc.vv vd, vs1, vd."""
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs2): CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (vd=vs2)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    # vd is both accumulator and vs2 source
    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    # vd serves as both accumulator and vs2
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC1}, {VREG_DST}")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_vd", format_data_line(vd_init, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Ternary accumulate: VV masked form (v0.t)
# ===================================================================

def emit_ternary_vv_masked(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vd_init: list[int],
    vs1: list[int],
    vs2: list[int],
    mask_bits: int,
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (masked)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_mask")
    tf.code(f"vlm.v v0, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC1}, {VREG_SRC2}, v0.t")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
    tf.data_align(sew)
    tf.data_label(f"{tag}_vd", format_data_line(vd_init, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Widening binary: VV masked form (v0.t)
# ===================================================================

def emit_widening_vv_masked(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    vs1: list[int],
    vd_init: list[int],
    mask_bits: int,
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    """Widening VV with v0.t mask. vd_init and expected are at 2*SEW."""
    dsew = 2 * sew
    nbytes = NUM_ELEMS * (dsew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew}\u2192e{dsew} {test_name} masked: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew}\u2192e{dsew} {test_name} masked: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew}\u2192e{dsew} {test_name} masked: CSR side-effect")
    tag = f"tc{cn_res}"
    src_nbytes = NUM_ELEMS * (sew // 8)

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (masked)")

    # Load vd_init at dest width
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
    tf.code(f"la t1, {tag}_vd")
    tf.code(f"vle{dsew}.v {VREG_DST}, (t1)")

    # Load sources at source width
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
    tf.code(f"la t1, {tag}_mask")
    tf.code(f"vlm.v v0, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")

    tf.code("SAVE_CSRS")
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0.t")

    # Check witness before vsetvli
    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {src_nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    # Store result at dest width
    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    wp = witness_pattern(sew)
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
    tf.data_align(dsew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_vd", format_data_line(vd_init, dsew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, dsew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Compare: VV masked form (v0.t — preserves inactive mask bits in vd)
# ===================================================================

def emit_compare_vv_masked(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    vs1: list[int],
    vd_init_mask: int,
    mask_bits: int,
    expected_mask: int,
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    """Compare VV with v0.t.  Inactive bits in vd keep vd_init_mask value."""
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: wrong mask result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} masked: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (masked)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")

    # Init vd (dest mask register) with known value
    tf.code(f"la t1, {tag}_vdinit")
    tf.code(f"vlm.v {VREG_DST}, (t1)")
    # Load active mask into v0
    tf.code(f"la t1, {tag}_mask")
    tf.code(f"vlm.v v0, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0.t")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vsm.v {VREG_DST}, (t1)")
    tf.code(f"lbu t2, 0(t1)")
    tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
    tf.code(f"li t3, {expected_mask}")
    tf.code("FAIL_IF_NE t2, t3")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
    tf.data_label(f"{tag}_vdinit", f"    .byte {vd_init_mask}")
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Generic binary: VV overlap form (vd == vs1)
# ===================================================================

def emit_binop_vv_overlap_vs1(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    vs2: list[int],
    vs1: list[int],
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    """Like emit_binop_vv but uses vd==vs1.  asm: {mnemonic} vd, vs2, vd."""
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs1): wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs1): witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} overlap(vd=vs1): CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (vd=vs1)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_s2")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    # Load vs1 into VREG_DST (overlap: vd == vs1)
    tf.code(f"la t1, {tag}_s1")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    # Execute: vd == vs1
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_DST}")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
    tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))


# ===================================================================
# Generic binary: VV self-overlap (vd == vs2 == vs1)
# ===================================================================

def emit_binop_vv_self_overlap(
    tf: TestFile,
    mnemonic: str,
    sew: int,
    test_name: str,
    src: list[int],
    expected: list[int],
    *,
    csr_check: str = "CHECK_CSRS_UNCHANGED",
) -> None:
    """VV with vd==vs2==vs1 (all same register).  asm: {mnemonic} vd, vd, vd."""
    nbytes = NUM_ELEMS * (sew // 8)
    cn_res = tf.next_check(f"{mnemonic} e{sew} {test_name} self-overlap: wrong result")
    cn_wit = tf.next_check(f"{mnemonic} e{sew} {test_name} self-overlap: witness changed")
    cn_csr = tf.next_check(f"{mnemonic} e{sew} {test_name} self-overlap: CSR side-effect")
    tag = f"tc{cn_res}"

    tf.blank()
    tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {test_name} (vd=vs2=vs1)")

    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")

    tf.code(f"la t1, {tag}_src")
    tf.code(f"vle{sew}.v {VREG_DST}, (t1)")

    tf.code(f"la t1, {tag}_w")
    tf.code(f"vle{sew}.v {VREG_WITNESS}, (t1)")
    tf.code("SAVE_CSRS")

    # Execute: all three operands are the same register
    tf.code(f"{mnemonic} {VREG_DST}, {VREG_DST}, {VREG_DST}")

    tf.code(f"SET_TEST_NUM {cn_res}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_wit}")
    tf.code(f"la t1, witness_buf")
    tf.code(f"vse{sew}.v {VREG_WITNESS}, (t1)")
    tf.code(f"CHECK_MEM witness_buf, {tag}_w, {nbytes}")

    tf.code(f"SET_TEST_NUM {cn_csr}")
    tf.code(csr_check)

    wp = witness_pattern(sew)
    tf.data_align(sew)
    tf.data_label(f"{tag}_src", format_data_line(src, sew))
    tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.data_label(f"{tag}_w", format_data_line(wp, sew))
