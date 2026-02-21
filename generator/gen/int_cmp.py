from __future__ import annotations

"""Generate tests for integer compare instructions:
vmseq, vmsne, vmsltu, vmslt, vmsleu, vmsle, vmsgtu, vmsgt.
"""

from pathlib import Path
from typing import Callable

from ..common import SEWS, U
from ..testfile import TestFile
from ..emit import emit_compare_vv, emit_compare_vv_masked
from ..vectors import compare_vv
from ..compute.integer import (
    mseq, msne, msltu, mslt, msleu, msle, msgtu, msgt,
)

_INSTRUCTIONS: list[tuple[str, list[str], Callable]] = [
    ("vmseq",  ["vv", "vx", "vi"], mseq),
    ("vmsne",  ["vv", "vx", "vi"], msne),
    ("vmsltu", ["vv", "vx"],       msltu),
    ("vmslt",  ["vv", "vx"],       mslt),
    ("vmsleu", ["vv", "vx", "vi"], msleu),
    ("vmsle",  ["vv", "vx", "vi"], msle),
    ("vmsgtu", ["vx", "vi"],       msgtu),
    ("vmsgt",  ["vx", "vi"],       msgt),
]


def _compute_mask(cfn: Callable, vs2: list[int], vs1: list[int], sew: int) -> int:
    """Compute expected mask from per-element compare function."""
    mask = 0
    for i, (a, b) in enumerate(zip(vs2, vs1)):
        if cfn(a, b, sew):
            mask |= 1 << i
    return mask


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []

    for base, forms, cfn in _INSTRUCTIONS:
        for form in forms:
            mnemonic = f"{base}.{form}"
            fpath = base_dir / "tests" / "int_cmp" / f"{base}_{form}.S"

            tf = TestFile(mnemonic, f"Integer compare {base} ({form} form)")

            if form == "vv":
                for sew in SEWS:
                    for name, s2, s1 in compare_vv(sew):
                        mask = _compute_mask(cfn, s2, s1, sew)
                        emit_compare_vv(tf, mnemonic, sew, name, s2, s1, mask)

                    # Masked compare at e32 only
                    if sew == 32:
                        s2_m = [5, 5, 1, 1]
                        s1_m = [5, 1, 5, 5]
                        vd_init_mask = 0b0101
                        mask_bits = 0b1010
                        expected_mask = 0
                        for i in range(4):
                            if (mask_bits >> i) & 1:
                                if cfn(s2_m[i], s1_m[i], sew):
                                    expected_mask |= 1 << i
                            else:
                                if (vd_init_mask >> i) & 1:
                                    expected_mask |= 1 << i
                        emit_compare_vv_masked(
                            tf, mnemonic, sew, "mixed",
                            s2_m, s1_m, vd_init_mask, mask_bits, expected_mask,
                        )

            elif form == "vx":
                from ..vectors import binop_vx
                for sew in SEWS:
                    for name, vec, sc in binop_vx(sew):
                        s1_expanded = [sc] * len(vec)
                        mask = _compute_mask(cfn, vec, s1_expanded, sew)
                        # Emit VX compare using inline code
                        cn_res = tf.next_check(f"{mnemonic} e{sew} {name}: wrong mask")
                        cn_csr = tf.next_check(f"{mnemonic} e{sew} {name}: CSR")
                        tag = f"tc{cn_res}"
                        from ..common import format_data_line, NUM_ELEMS
                        from ..emit import VREG_SRC2, VREG_DST

                        tf.blank()
                        tf.comment(f"Test {cn_res}: {mnemonic} SEW={sew} {name}")
                        tf.code(f"li t0, {NUM_ELEMS}")
                        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                        tf.code(f"la t1, {tag}_s2")
                        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                        from ..common import format_hex
                        tf.code(f"li a0, {format_hex(U(sc, sew), sew)}")
                        tf.code("SAVE_CSRS")
                        tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, a0")
                        tf.code(f"SET_TEST_NUM {cn_res}")
                        tf.code(f"la t1, result_buf")
                        tf.code(f"vsm.v {VREG_DST}, (t1)")
                        tf.code(f"lbu t2, 0(t1)")
                        tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
                        tf.code(f"li t3, {mask}")
                        tf.code("FAIL_IF_NE t2, t3")
                        tf.code(f"SET_TEST_NUM {cn_csr}")
                        tf.code("CHECK_CSRS_UNCHANGED")
                        tf.data_align(sew)
                        tf.data_label(f"{tag}_s2", format_data_line(vec, sew))

            elif form == "vi":
                from ..vectors import binop_vi
                for sew in SEWS:
                    for name, vec, imm in binop_vi(sew):
                        imm_val = U(imm, sew)
                        s1_expanded = [imm_val] * len(vec)
                        mask = _compute_mask(cfn, vec, s1_expanded, sew)
                        cn_res = tf.next_check(f"{mnemonic} e{sew} {name}: wrong mask")
                        cn_csr = tf.next_check(f"{mnemonic} e{sew} {name}: CSR")
                        tag = f"tc{cn_res}"
                        from ..common import format_data_line, NUM_ELEMS
                        from ..emit import VREG_SRC2, VREG_DST

                        tf.blank()
                        tf.comment(f"Test {cn_res}: {mnemonic} SEW={sew} {name}")
                        tf.code(f"li t0, {NUM_ELEMS}")
                        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                        tf.code(f"la t1, {tag}_s2")
                        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                        tf.code("SAVE_CSRS")
                        tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {imm}")
                        tf.code(f"SET_TEST_NUM {cn_res}")
                        tf.code(f"la t1, result_buf")
                        tf.code(f"vsm.v {VREG_DST}, (t1)")
                        tf.code(f"lbu t2, 0(t1)")
                        tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
                        tf.code(f"li t3, {mask}")
                        tf.code("FAIL_IF_NE t2, t3")
                        tf.code(f"SET_TEST_NUM {cn_csr}")
                        tf.code("CHECK_CSRS_UNCHANGED")
                        tf.data_align(sew)
                        tf.data_label(f"{tag}_s2", format_data_line(vec, sew))

            tf.write(fpath)
            generated.append(str(fpath))

    return generated
