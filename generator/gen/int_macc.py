from __future__ import annotations

"""Generate tests for integer multiply-accumulate instructions:
vmacc, vnmsac, vmadd, vnmsub (vv/vx forms),
vwmaccu, vwmacc, vwmaccsu, vwmaccus (vv/vx forms).
"""

from pathlib import Path
from typing import Callable

from ..common import SEWS, WIDENING_SEWS, M, U, format_data_line, format_hex, NUM_ELEMS
from ..testfile import TestFile
from ..emit import (
    emit_ternary_vv, emit_ternary_vv_overlap, emit_ternary_vv_masked,
    VREG_DST, VREG_SRC1, VREG_SRC2,
)
from ..vectors import macc_vv
from ..compute.integer import (
    macc, nmsac, madd, nmsub,
    wmaccu, wmacc, wmaccsu, wmaccus,
)

_SINGLE_WIDTH: list[tuple[str, Callable]] = [
    ("vmacc.vv",  macc),
    ("vnmsac.vv", nmsac),
    ("vmadd.vv",  madd),
    ("vnmsub.vv", nmsub),
]


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "int_macc"

    # --- Single-width VV ---
    for mnemonic, cfn in _SINGLE_WIDTH:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Multiply-accumulate {mnemonic}")
        for sew in SEWS:
            for name, vd_init, vs1, vs2 in macc_vv(sew):
                exp = [cfn(d, s1, s2, sew)
                       for d, s1, s2 in zip(vd_init, vs1, vs2)]
                emit_ternary_vv(tf, mnemonic, sew, name,
                                vd_init, vs1, vs2, exp)

            # Overlap (vd==vs2) + Masked at e32 only
            if sew == 32:
                vd_ov = [10, 20, 30, 40]
                vs1_ov = [1, 2, 3, 4]
                # vd serves as both accumulator and vs2
                exp_ov = [cfn(d, s1, d, sew)
                          for d, s1 in zip(vd_ov, vs1_ov)]
                emit_ternary_vv_overlap(
                    tf, mnemonic, sew, "basic", vd_ov, vs1_ov, exp_ov,
                )
                # Masked test
                vd_init_m = [0xDEAD] * 4
                vs2_m = [10, 20, 30, 40]
                mask_bits = 0b1010
                exp_full = [cfn(d, s1, s2, sew)
                            for d, s1, s2 in zip(vd_init_m, vs1_ov, vs2_m)]
                exp_masked = [
                    exp_full[i] if (mask_bits >> i) & 1 else vd_init_m[i]
                    for i in range(4)
                ]
                emit_ternary_vv_masked(
                    tf, mnemonic, sew, "basic",
                    vd_init_m, vs1_ov, vs2_m, mask_bits, exp_masked,
                )
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Single-width VX ---
    for base_mn, cfn in _SINGLE_WIDTH:
        mn_vx = base_mn.replace(".vv", ".vx")
        fname = mn_vx.replace(".", "_") + ".S"
        tf = TestFile(mn_vx, f"Multiply-accumulate {mn_vx}")
        for sew in SEWS:
            m = M(sew)
            test_cases = [
                ("basic", [0, 0, 0, 0], 2, [10, 20, 30, 40]),
                ("accum", [100, 200, 300, 400], 1, [1, 2, 3, 4]),
            ]
            for name, vd_init, scalar, vs2 in test_cases:
                exp = [cfn(d, scalar, s2, sew)
                       for d, s2 in zip(vd_init, vs2)]
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mn_vx} e{sew} {name}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.comment(f"Test {cn_res}: {mn_vx} SEW={sew} {name}")
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code(f"la t1, {tag}_vd")
                tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
                tf.code(f"la t1, {tag}_s2")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                tf.code(f"li a0, {format_hex(U(scalar, sew), sew)}")
                tf.code("SAVE_CSRS")
                tf.code(f"{mn_vx} {VREG_DST}, a0, {VREG_SRC2}")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(sew)
                tf.data_label(f"{tag}_vd", format_data_line(vd_init, sew))
                tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Widening multiply-accumulate VV ---
    _WMACC = [
        ("vwmaccu.vv", wmaccu),
        ("vwmacc.vv",  wmacc),
        ("vwmaccsu.vv", wmaccsu),
    ]
    for mnemonic, cfn in _WMACC:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Widening macc {mnemonic}")
        for sew in WIDENING_SEWS:
            dsew = 2 * sew
            vd_init = [0] * NUM_ELEMS
            vs1 = [1, 2, 3, 4]
            vs2 = [10, 20, 30, 40]
            exp = [cfn(d, s1, s2, sew) for d, s1, s2 in zip(vd_init, vs1, vs2)]
            nbytes = NUM_ELEMS * (dsew // 8)
            cn_res = tf.next_check(f"{mnemonic} e{sew} basic: result")
            cn_csr = tf.next_check(f"{mnemonic} e{sew} basic: CSR side-effect")
            tag = f"tc{cn_res}"

            tf.blank()
            tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew}")
            # Load vd at dest width
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, {tag}_vd")
            tf.code(f"vle{dsew}.v {VREG_DST}, (t1)")
            # Load sources at source width
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC1}, {VREG_SRC2}")
            # Check CSRs before vsetvli changes vl/vtype
            tf.code(f"SET_TEST_NUM {cn_csr}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.code(f"SET_TEST_NUM {cn_res}")
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

            tf.data_align(dsew)
            tf.data_label(f"{tag}_vd", format_data_line(vd_init, dsew))
            tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
            tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Widening multiply-accumulate VX ---
    _WMACC_VX = [
        ("vwmaccu.vx", wmaccu),
        ("vwmacc.vx",  wmacc),
        ("vwmaccsu.vx", wmaccsu),
        ("vwmaccus.vx", wmaccus),
    ]
    for mnemonic, cfn in _WMACC_VX:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Widening macc {mnemonic}")
        for sew in WIDENING_SEWS:
            dsew = 2 * sew
            vd_init = [0] * NUM_ELEMS
            vs2 = [10, 20, 30, 40]
            scalar = 2
            exp = [cfn(d, scalar, s2, sew) for d, s2 in zip(vd_init, vs2)]
            nbytes = NUM_ELEMS * (dsew // 8)
            cn_res = tf.next_check(f"{mnemonic} e{sew} basic: result")
            cn_csr = tf.next_check(f"{mnemonic} e{sew} basic: CSR side-effect")
            tag = f"tc{cn_res}"

            tf.blank()
            tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew}")
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, {tag}_vd")
            tf.code(f"vle{dsew}.v {VREG_DST}, (t1)")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"li a0, {format_hex(U(scalar, sew), sew)}")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, a0, {VREG_SRC2}")
            # Check CSRs before vsetvli changes vl/vtype
            tf.code(f"SET_TEST_NUM {cn_csr}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.code(f"SET_TEST_NUM {cn_res}")
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")

            tf.data_align(dsew)
            tf.data_label(f"{tag}_vd", format_data_line(vd_init, dsew))
            tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    return generated
