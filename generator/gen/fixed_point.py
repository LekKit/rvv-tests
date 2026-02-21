from __future__ import annotations

"""Generate tests for fixed-point instructions:
vsaddu, vsadd, vssubu, vssub, vaaddu, vaadd, vasubu, vasub,
vsmul, vssrl, vssra, vnclipu, vnclip.
"""

from pathlib import Path

from ..common import SEWS, WIDENING_SEWS, M, U, format_data_line, format_hex, NUM_ELEMS
from ..testfile import TestFile
from ..emit import (
    emit_binop_vv, emit_binop_vv_overlap, emit_binop_vv_masked,
    VREG_SRC2, VREG_SRC1, VREG_DST,
)
from ..vectors import binop_vv, binop_vx
from ..compute.fixed_point import (
    saddu, sadd, ssubu, ssub,
    aaddu, aadd, asubu, asub,
    smul, ssrl, ssra,
    nclipu, nclip,
)


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "fixed_point"

    # --- Saturating add/sub VV ---
    _SAT_OPS = [
        ("vsaddu.vv", saddu), ("vsadd.vv", sadd),
        ("vssubu.vv", ssubu), ("vssub.vv", ssub),
    ]
    for mnemonic, cfn in _SAT_OPS:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Saturating {mnemonic}")
        for sew in SEWS:
            for name, s2, s1 in binop_vv(sew):
                results = [cfn(a, b, sew) for a, b in zip(s2, s1)]
                exp = [r[0] for r in results]
                emit_binop_vv(tf, mnemonic, sew, name, s2, s1, exp,
                              csr_check="CHECK_CSRS_UNCHANGED_FIXEDPOINT")

            # Overlap + Masked at e32 only
            if sew == 32:
                s2_ov, s1_ov = [1, 2, 3, 4], [10, 20, 30, 40]
                results_ov = [cfn(a, b, sew) for a, b in zip(s2_ov, s1_ov)]
                exp_ov = [r[0] for r in results_ov]
                emit_binop_vv_overlap(
                    tf, mnemonic, sew, "basic", s2_ov, s1_ov, exp_ov,
                    csr_check="CHECK_CSRS_UNCHANGED_FIXEDPOINT",
                )
                vd_init = [0xDEAD] * 4
                mask_bits = 0b1010
                exp_masked = [
                    exp_ov[i] if (mask_bits >> i) & 1 else vd_init[i]
                    for i in range(4)
                ]
                emit_binop_vv_masked(
                    tf, mnemonic, sew, "basic",
                    s2_ov, s1_ov, vd_init, mask_bits, exp_masked,
                    csr_check="CHECK_CSRS_UNCHANGED_FIXEDPOINT",
                )
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Averaging add/sub VV ---
    _AVG_OPS = [
        ("vaaddu.vv", aaddu), ("vaadd.vv", aadd),
        ("vasubu.vv", asubu), ("vasub.vv", asub),
    ]
    for mnemonic, cfn in _AVG_OPS:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Averaging {mnemonic}")
        for sew in SEWS:
            for name, s2, s1 in binop_vv(sew):
                # Test with vxrm=0 (round-to-nearest-up)
                exp = [cfn(a, b, sew, vxrm=0) for a, b in zip(s2, s1)]
                # Need to set vxrm before the instruction
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mnemonic} e{sew} {name}: result")
                cn_csr = tf.next_check(f"{mnemonic} e{sew} {name}: CSR")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.comment(f"Test {cn_res}: {mnemonic} SEW={sew} {name} vxrm=0")
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code("csrwi vxrm, 0")
                tf.code(f"la t1, {tag}_s2")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                tf.code(f"la t1, {tag}_s1")
                tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code(f"SET_TEST_NUM {cn_csr}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(sew)
                tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
                tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Scaling shift VV ---
    _SSHIFT_OPS = [("vssrl.vv", ssrl), ("vssra.vv", ssra)]
    for mnemonic, cfn in _SSHIFT_OPS:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Scaling shift {mnemonic}")
        for sew in SEWS:
            m = M(sew)
            test_cases = [
                ("shift0", [m, 1, 0, m >> 1], [0, 0, 0, 0]),
                ("shift1", [m, 0x55 & m, 0xAA & m, 1], [1, 1, 1, 1]),
            ]
            for name, s2, s1 in test_cases:
                exp = [cfn(a, b, sew, vxrm=0) for a, b in zip(s2, s1)]
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mnemonic} e{sew} {name}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code("csrwi vxrm, 0")
                tf.code(f"la t1, {tag}_s2")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                tf.code(f"la t1, {tag}_s1")
                tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(sew)
                tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
                tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Saturating add/sub VX ---
    _SAT_VX = [
        ("vsaddu.vx", saddu), ("vsadd.vx", sadd),
        ("vssubu.vx", ssubu), ("vssub.vx", ssub),
    ]
    for mnemonic, cfn in _SAT_VX:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Saturating {mnemonic}")
        for sew in SEWS:
            for name, s2, sc in binop_vx(sew):
                results = [cfn(a, sc, sew) for a in s2]
                exp = [r[0] for r in results]
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mnemonic} e{sew} {name}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code(f"la t1, {tag}_s2")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                tf.code(f"li a0, {format_hex(sc, sew)}")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, a0")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED_FIXEDPOINT")

                tf.data_align(sew)
                tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Saturating add VI (only saddu and sadd have VI forms) ---
    _SAT_VI = [("vsaddu.vi", saddu), ("vsadd.vi", sadd)]
    for mnemonic, cfn in _SAT_VI:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Saturating {mnemonic}")
        for sew in SEWS:
            m = M(sew)
            for imm in [0, 1, 15, -1, -16]:
                s2 = [0, 1, m, m >> 1]
                results = [cfn(a, imm, sew) for a in s2]
                exp = [r[0] for r in results]
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mnemonic} e{sew} imm={imm}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code(f"la t1, {tag}_s2")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {imm}")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED_FIXEDPOINT")

                tf.data_align(sew)
                tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Averaging add/sub VX ---
    _AVG_VX = [
        ("vaaddu.vx", aaddu), ("vaadd.vx", aadd),
        ("vasubu.vx", asubu), ("vasub.vx", asub),
    ]
    for mnemonic, cfn in _AVG_VX:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Averaging {mnemonic}")
        for sew in SEWS:
            for name, s2, sc in binop_vx(sew):
                exp = [cfn(a, sc, sew, vxrm=0) for a in s2]
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mnemonic} e{sew} {name}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code("csrwi vxrm, 0")
                tf.code(f"la t1, {tag}_s2")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                tf.code(f"li a0, {format_hex(sc, sew)}")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, a0")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(sew)
                tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- vsmul.vv / vsmul.vx ---
    for form in ["vv", "vx"]:
        mnemonic = f"vsmul.{form}"
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Fractional multiply {mnemonic}")
        for sew in SEWS:
            m = M(sew)
            h = 1 << (sew - 1)
            if form == "vv":
                cases: list[tuple[str, list[int], list[int] | int]] = [
                    ("basic", [h, h + 1, m, 1], [h, h + 1, m, 1]),
                    ("mixed", [0, 1, h - 1, m], [m, h, 1, 0]),
                ]
            else:
                cases = [
                    ("basic", [h, h + 1, m, 1], h),
                    ("one",   [0, 1, h - 1, m], 1),
                ]
            for name, s2, s1_or_sc in cases:
                if form == "vv":
                    results = [smul(a, b, sew, vxrm=0) for a, b in zip(s2, s1_or_sc)]
                else:
                    results = [smul(a, s1_or_sc, sew, vxrm=0) for a in s2]
                exp = [r[0] for r in results]
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mnemonic} e{sew} {name}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code("csrwi vxrm, 0")
                tf.code(f"la t1, {tag}_s2")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                if form == "vv":
                    tf.code(f"la t1, {tag}_s1")
                    tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
                else:
                    tf.code(f"li a0, {format_hex(s1_or_sc, sew)}")
                tf.code("SAVE_CSRS")
                if form == "vv":
                    tf.code(f"vsmul.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
                else:
                    tf.code(f"vsmul.vx {VREG_DST}, {VREG_SRC2}, a0")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED_FIXEDPOINT")

                tf.data_align(sew)
                tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
                if form == "vv":
                    tf.data_label(f"{tag}_s1", format_data_line(s1_or_sc, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Scaling shift VX ---
    _SSHIFT_VX = [("vssrl.vx", ssrl), ("vssra.vx", ssra)]
    for mnemonic, cfn in _SSHIFT_VX:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Scaling shift {mnemonic}")
        for sew in SEWS:
            m = M(sew)
            for shift_amt in [0, 1, sew - 1]:
                s2 = [m, 0x55 & m, 0xAA & m, 1]
                exp = [cfn(a, shift_amt, sew, vxrm=0) for a in s2]
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mnemonic} e{sew} sh={shift_amt}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code("csrwi vxrm, 0")
                tf.code(f"la t1, {tag}_s2")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                tf.code(f"li a0, {shift_amt}")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, a0")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(sew)
                tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Scaling shift VI ---
    _SSHIFT_VI = [("vssrl.vi", ssrl), ("vssra.vi", ssra)]
    for mnemonic, cfn in _SSHIFT_VI:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Scaling shift {mnemonic}")
        for sew in SEWS:
            m = M(sew)
            for imm in [0, 1, min(7, sew - 1), min(15, sew - 1)]:
                s2 = [m, 0x55 & m, 0xAA & m, 1]
                exp = [cfn(a, imm, sew, vxrm=0) for a in s2]
                nbytes = NUM_ELEMS * (sew // 8)
                cn_res = tf.next_check(f"{mnemonic} e{sew} imm={imm}: result")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code("csrwi vxrm, 0")
                tf.code(f"la t1, {tag}_s2")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {imm}")
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(sew)
                tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Narrowing clip: vnclipu / vnclip  (WV/WX/WI) ---
    _NCLIP = [("vnclipu", nclipu), ("vnclip", nclip)]
    for base, cfn in _NCLIP:
        for form in ["wv", "wx", "wi"]:
            mnemonic = f"{base}.{form}"
            fname = mnemonic.replace(".", "_") + ".S"
            tf = TestFile(mnemonic, f"Narrowing clip {mnemonic}")
            for sew in WIDENING_SEWS:
                dsew = 2 * sew
                dm = M(dsew)
                m = M(sew)
                if form == "wv":
                    cases_n = [
                        ("basic", [0x100, 0x200, dm, 0], [0, 0, 0, 0]),
                        ("shift", [dm, dm >> 1, 1, 0], [sew, sew // 2, 0, 1]),
                    ]
                elif form == "wx":
                    cases_n = [
                        ("sh0", [0x100, 0x200, dm, 0], 0),
                        ("sh1", [dm, dm >> 1, 1, 0], 1),
                    ]
                else:  # wi
                    cases_n = [
                        ("imm0", [0x100, 0x200, dm, 0], 0),
                        ("imm1", [dm, dm >> 1, 1, 0], 1),
                    ]
                for name, s2_wide, s1_or_sc in cases_n:
                    if form == "wv":
                        results = [cfn(a, b, sew, vxrm=0) for a, b in zip(s2_wide, s1_or_sc)]
                    else:
                        results = [cfn(a, s1_or_sc, sew, vxrm=0) for a in s2_wide]
                    exp = [r[0] for r in results]
                    nbytes = NUM_ELEMS * (sew // 8)
                    cn_res = tf.next_check(f"{mnemonic} e{sew} {name}: result")
                    tag = f"tc{cn_res}"

                    tf.blank()
                    tf.code(f"li t0, {NUM_ELEMS}")
                    tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
                    tf.code("csrwi vxrm, 0")
                    tf.code(f"la t1, {tag}_s2w")
                    tf.code(f"vle{dsew}.v {VREG_SRC2}, (t1)")
                    if form == "wv":
                        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                        tf.code(f"la t1, {tag}_s1")
                        tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
                    elif form == "wx":
                        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                        tf.code(f"li a0, {s1_or_sc}")
                    else:
                        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                    tf.code("SAVE_CSRS")
                    if form == "wv":
                        tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
                    elif form == "wx":
                        tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, a0")
                    else:
                        tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {s1_or_sc}")
                    tf.code(f"SET_TEST_NUM {cn_res}")
                    tf.code(f"la t1, result_buf")
                    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                    tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                    tf.code("CHECK_CSRS_UNCHANGED_FIXEDPOINT")

                    tf.data_align(dsew)
                    tf.data_label(f"{tag}_s2w", format_data_line(s2_wide, dsew))
                    if form == "wv":
                        tf.data_label(f"{tag}_s1", format_data_line(s1_or_sc, sew))
                    tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
            fpath = out / fname
            tf.write(fpath)
            generated.append(str(fpath))

    return generated
