from __future__ import annotations

"""Generate tests for widening/narrowing integer instructions:
vwaddu, vwadd, vwsubu, vwsub (.vv/.vx/.wv/.wx),
vwmulu, vwmul, vwmulsu (.vv/.vx),
vnsrl, vnsra (.wv/.wx/.wi).
"""

from pathlib import Path
from typing import Callable

from ..common import WIDENING_SEWS, SEWS, M, U
from ..testfile import TestFile
from ..emit import emit_widening_vv, emit_widening_vv_masked, emit_binop_vv
from ..vectors import widening_vv, widening_wv, narrowing_wv
from ..compute.integer import (
    waddu, wadd, wsubu, wsub,
    waddu_w, wadd_w, wsubu_w, wsub_w,
    wmulu, wmul, wmulsu,
    nsrl, nsra,
)

_WIDENING_VV: list[tuple[str, Callable]] = [
    ("vwaddu.vv", waddu),
    ("vwadd.vv",  wadd),
    ("vwsubu.vv", wsubu),
    ("vwsub.vv",  wsub),
    ("vwmulu.vv", wmulu),
    ("vwmul.vv",  wmul),
    ("vwmulsu.vv", wmulsu),
]

_WIDENING_WV: list[tuple[str, Callable]] = [
    ("vwaddu.wv", waddu_w),
    ("vwadd.wv",  wadd_w),
    ("vwsubu.wv", wsubu_w),
    ("vwsub.wv",  wsub_w),
]

_NARROWING_WV: list[tuple[str, Callable]] = [
    ("vnsrl.wv", nsrl),
    ("vnsra.wv", nsra),
]


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "int_widening"

    # --- Widening VV ---
    for mnemonic, cfn in _WIDENING_VV:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Widening {mnemonic}")
        for sew in WIDENING_SEWS:
            dsew = 2 * sew
            for name, s2, s1 in widening_vv(sew):
                exp = [cfn(a, b, sew) for a, b in zip(s2, s1)]
                emit_widening_vv(tf, mnemonic, sew, name, s2, s1, exp)

            # Masked widening at sew == 16 only
            if sew == 16:
                s2_m = [1, 2, 3, 4]
                s1_m = [5, 6, 7, 8]
                vd_init_m = [0xDEAD] * 4  # at dsew (32-bit)
                mask_bits = 0b1010
                exp_full = [cfn(a, b, sew) for a, b in zip(s2_m, s1_m)]
                exp_masked = [
                    exp_full[i] if (mask_bits >> i) & 1 else vd_init_m[i]
                    for i in range(4)
                ]
                emit_widening_vv_masked(
                    tf, mnemonic, sew, "basic",
                    s2_m, s1_m, vd_init_m, mask_bits, exp_masked,
                )
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Widening VX ---
    for base_mn, cfn in _WIDENING_VV:
        mn_vx = base_mn.replace(".vv", ".vx")
        fname = mn_vx.replace(".", "_") + ".S"
        tf = TestFile(mn_vx, f"Widening {mn_vx}")
        for sew in WIDENING_SEWS:
            dsew = 2 * sew
            from ..common import format_data_line, format_hex, NUM_ELEMS
            from ..emit import VREG_SRC2, VREG_DST
            m = M(sew)
            scalars = [1, m, 0, 1 << (sew - 1)]
            for sc in scalars:
                name = f"sc_{format_hex(sc, sew)}"
                vec = [1, 2, m, m >> 1]
                exp = [cfn(a, sc, sew) for a in vec]
                nbytes = NUM_ELEMS * (dsew // 8)
                cn_res = tf.next_check(f"{mn_vx} e{sew} {name}: result")
                cn_csr = tf.next_check(f"{mn_vx} e{sew} {name}: CSR side-effect")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.comment(f"Test {cn_res}-{cn_csr}: {mn_vx} SEW={sew} {name}")
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code(f"la t1, {tag}_s2")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                tf.code(f"li a0, {format_hex(sc, sew)}")
                tf.code("SAVE_CSRS")
                tf.code(f"{mn_vx} {VREG_DST}, {VREG_SRC2}, a0")
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
                tf.data_label(f"{tag}_s2", format_data_line(vec, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Widening .W V ---
    for mnemonic, cfn in _WIDENING_WV:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Widening-W {mnemonic}")
        for sew in WIDENING_SEWS:
            dsew = 2 * sew
            for name, s2_wide, s1 in widening_wv(sew):
                exp = [cfn(a, b, sew) for a, b in zip(s2_wide, s1)]
                from ..common import format_data_line, NUM_ELEMS
                from ..emit import VREG_DST, VREG_SRC2, VREG_SRC1
                nbytes = NUM_ELEMS * (dsew // 8)
                cn_res = tf.next_check(f"{mnemonic} e{sew} {name}: result")
                cn_csr = tf.next_check(f"{mnemonic} e{sew} {name}: CSR side-effect")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.comment(f"Test {cn_res}-{cn_csr}: {mnemonic} SEW={sew} {name}")
                # Load wide source (2*SEW)
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
                tf.code(f"la t1, {tag}_s2w")
                tf.code(f"vle{dsew}.v {VREG_SRC2}, (t1)")
                # Load narrow source (SEW)
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code(f"la t1, {tag}_s1")
                tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
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
                tf.data_label(f"{tag}_s2w", format_data_line(s2_wide, dsew))
                tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Narrowing WV ---
    for mnemonic, cfn in _NARROWING_WV:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Narrowing {mnemonic}")
        for sew in WIDENING_SEWS:
            dsew = 2 * sew
            for name, s2_wide, s1 in narrowing_wv(sew):
                exp = [cfn(a, b, sew) for a, b in zip(s2_wide, s1)]
                from ..common import format_data_line, NUM_ELEMS
                from ..emit import VREG_DST, VREG_SRC2, VREG_SRC1
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mnemonic} e{sew} {name}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.comment(f"Test {cn_res}: {mnemonic} SEW={sew} {name}")
                # Load wide source at 2*SEW
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
                tf.code(f"la t1, {tag}_s2w")
                tf.code(f"vle{dsew}.v {VREG_SRC2}, (t1)")
                # Load shift amount at SEW
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code(f"la t1, {tag}_s1")
                tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(dsew)
                tf.data_label(f"{tag}_s2w", format_data_line(s2_wide, dsew))
                tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Narrowing WX ---
    for mnemonic, cfn in _NARROWING_WV:
        mn_wx = mnemonic.replace(".wv", ".wx")
        fname = mn_wx.replace(".", "_") + ".S"
        tf = TestFile(mn_wx, f"Narrowing {mn_wx}")
        for sew in WIDENING_SEWS:
            dsew = 2 * sew
            from ..common import format_data_line, NUM_ELEMS
            from ..emit import VREG_DST, VREG_SRC2
            dm = M(dsew)
            for shift_amt in [0, 1, sew]:
                s2_wide = [dm, dm >> 1, 1, 0]
                exp = [cfn(a, shift_amt, sew) for a in s2_wide]
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mn_wx} e{sew} sh={shift_amt}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
                tf.code(f"la t1, {tag}_s2w")
                tf.code(f"vle{dsew}.v {VREG_SRC2}, (t1)")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code(f"li a0, {shift_amt}")
                tf.code("SAVE_CSRS")
                tf.code(f"{mn_wx} {VREG_DST}, {VREG_SRC2}, a0")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(dsew)
                tf.data_label(f"{tag}_s2w", format_data_line(s2_wide, dsew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Narrowing WI ---
    for mnemonic, cfn in _NARROWING_WV:
        mn_wi = mnemonic.replace(".wv", ".wi")
        fname = mn_wi.replace(".", "_") + ".S"
        tf = TestFile(mn_wi, f"Narrowing {mn_wi}")
        for sew in WIDENING_SEWS:
            dsew = 2 * sew
            from ..common import format_data_line, NUM_ELEMS
            from ..emit import VREG_DST, VREG_SRC2
            dm = M(dsew)
            for imm in [0, 1, min(15, sew)]:
                s2_wide = [dm, dm >> 1, 1, 0]
                exp = [cfn(a, imm, sew) for a in s2_wide]
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mn_wi} e{sew} imm={imm}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
                tf.code(f"la t1, {tag}_s2w")
                tf.code(f"vle{dsew}.v {VREG_SRC2}, (t1)")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code("SAVE_CSRS")
                tf.code(f"{mn_wi} {VREG_DST}, {VREG_SRC2}, {imm}")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(dsew)
                tf.data_label(f"{tag}_s2w", format_data_line(s2_wide, dsew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Widening .WX (wide + scalar) ---
    for mnemonic, cfn in _WIDENING_WV:
        mn_wx = mnemonic.replace(".wv", ".wx")
        fname = mn_wx.replace(".", "_") + ".S"
        tf = TestFile(mn_wx, f"Widening-W {mn_wx}")
        for sew in WIDENING_SEWS:
            dsew = 2 * sew
            from ..common import format_data_line, format_hex, NUM_ELEMS
            from ..emit import VREG_DST, VREG_SRC2
            dm = M(dsew)
            m = M(sew)
            for sc in [1, m, 0]:
                s2_wide = [1, 2, dm >> 1, dm]
                exp = [cfn(a, sc, sew) for a in s2_wide]
                nbytes = NUM_ELEMS * (dsew // 8)
                name = f"sc_{format_hex(U(sc, sew), sew)}"
                cn_res = tf.next_check(f"{mn_wx} e{sew} {name}: result")
                cn_csr = tf.next_check(f"{mn_wx} e{sew} {name}: CSR side-effect")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
                tf.code(f"la t1, {tag}_s2w")
                tf.code(f"vle{dsew}.v {VREG_SRC2}, (t1)")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code(f"li a0, {format_hex(sc, sew)}")
                tf.code("SAVE_CSRS")
                tf.code(f"{mn_wx} {VREG_DST}, {VREG_SRC2}, a0")
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
                tf.data_label(f"{tag}_s2w", format_data_line(s2_wide, dsew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    return generated
