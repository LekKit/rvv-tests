from __future__ import annotations

"""Generate tests for floating-point arithmetic instructions:
vfadd, vfsub, vfrsub, vfmul, vfdiv, vfrdiv,
vfmin, vfmax, vfsgnj, vfsgnjn, vfsgnjx,
vfmacc, vfnmacc, vfmsac, vfnmsac, vfmadd, vfnmadd, vfmsub, vfnmsub,
vfsqrt, vfclass,
vmfeq, vmfne, vmflt, vmfle, vmfgt, vmfge,
vfcvt.*, vfmerge, vfmv.
"""

from pathlib import Path
from typing import Callable

from ..common import FP_SEWS, NUM_ELEMS, format_data_line
from ..testfile import TestFile
from ..emit import (
    emit_fp_binop_vv, emit_fp_binop_vf, emit_compare_vv, emit_compare_vf,
    emit_ternary_vv, emit_ternary_vf,
    emit_binop_vv_overlap, emit_binop_vv_masked,
    emit_fp_binop_vf_overlap, emit_fp_binop_vf_masked,
    emit_ternary_vv_overlap, emit_ternary_vv_masked,
    emit_compare_vv_masked,
    VREG_DST, VREG_SRC2, VREG_SRC1,
)
from ..vectors import fp_binop_vv, fp_binop_vf, fp_fma_vv, fp_fma_vf
from ..compute.floating_point import (
    fadd, fsub, frsub, fmul, fdiv, frdiv,
    fmin, fmax, fsgnj, fsgnjn, fsgnjx,
    mfeq, mfne, mflt, mfle, mfgt, mfge,
    fclass, fsqrt,
    fmacc, fnmacc, fmsac, fnmsac, fmadd, fnmadd, fmsub, fnmsub,
    fcvt_xu_f, fcvt_x_f, fcvt_rtz_xu_f, fcvt_rtz_x_f, fcvt_f_xu, fcvt_f_x,
)

_FP_BINOP_VV: list[tuple[str, Callable, str]] = [
    ("vfadd.vv",   fadd,   "float_arith"),
    ("vfsub.vv",   fsub,   "float_arith"),
    ("vfmul.vv",   fmul,   "float_arith"),
    ("vfdiv.vv",   fdiv,   "float_arith"),
    ("vfmin.vv",   fmin,   "float_minmax"),
    ("vfmax.vv",   fmax,   "float_minmax"),
    ("vfsgnj.vv",  fsgnj,  "float_sgnj"),
    ("vfsgnjn.vv", fsgnjn, "float_sgnj"),
    ("vfsgnjx.vv", fsgnjx, "float_sgnj"),
]

_FP_BINOP_VF: list[tuple[str, Callable, str]] = [
    ("vfadd.vf",   fadd,   "float_arith"),
    ("vfsub.vf",   fsub,   "float_arith"),
    ("vfrsub.vf",  frsub,  "float_arith"),
    ("vfmul.vf",   fmul,   "float_arith"),
    ("vfdiv.vf",   fdiv,   "float_arith"),
    ("vfrdiv.vf",  frdiv,  "float_arith"),
    ("vfmin.vf",   fmin,   "float_minmax"),
    ("vfmax.vf",   fmax,   "float_minmax"),
    ("vfsgnj.vf",  fsgnj,  "float_sgnj"),
    ("vfsgnjn.vf", fsgnjn, "float_sgnj"),
    ("vfsgnjx.vf", fsgnjx, "float_sgnj"),
]

_FP_FMA_VV: list[tuple[str, Callable]] = [
    ("vfmacc.vv",  fmacc),
    ("vfnmacc.vv", fnmacc),
    ("vfmsac.vv",  fmsac),
    ("vfnmsac.vv", fnmsac),
    ("vfmadd.vv",  fmadd),
    ("vfnmadd.vv", fnmadd),
    ("vfmsub.vv",  fmsub),
    ("vfnmsub.vv", fnmsub),
]

_FP_FMA_VF: list[tuple[str, Callable]] = [
    ("vfmacc.vf",  fmacc),
    ("vfnmacc.vf", fnmacc),
    ("vfmsac.vf",  fmsac),
    ("vfnmsac.vf", fnmsac),
    ("vfmadd.vf",  fmadd),
    ("vfnmadd.vf", fnmadd),
    ("vfmsub.vf",  fmsub),
    ("vfnmsub.vf", fnmsub),
]

_FP_CMP_VV: list[tuple[str, Callable]] = [
    ("vmfeq.vv", mfeq), ("vmfne.vv", mfne),
    ("vmflt.vv", mflt), ("vmfle.vv", mfle),
]

_FP_CMP_VF: list[tuple[str, Callable]] = [
    ("vmfeq.vf", mfeq), ("vmfne.vf", mfne),
    ("vmflt.vf", mflt), ("vmfle.vf", mfle),
    ("vmfgt.vf", mfgt), ("vmfge.vf", mfge),
]


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []

    # --- FP Binary VV ---
    from ..common import f32_to_bits as _b32
    for mnemonic, cfn, subdir in _FP_BINOP_VV:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / subdir / fname
        tf = TestFile(mnemonic, f"FP {mnemonic}")
        for sew in FP_SEWS:
            for name, s2, s1 in fp_binop_vv(sew):
                exp = [cfn(a, b, sew) for a, b in zip(s2, s1)]
                emit_fp_binop_vv(tf, mnemonic, sew, name, s2, s1, exp)

            # Overlap + Masked at e32 only
            if sew == 32:
                s2_ov = [_b32(1.0), _b32(2.0), _b32(3.0), _b32(4.0)]
                s1_ov = [_b32(0.5), _b32(0.25), _b32(0.125), _b32(0.0625)]
                exp_ov = [cfn(a, b, sew) for a, b in zip(s2_ov, s1_ov)]
                emit_binop_vv_overlap(
                    tf, mnemonic, sew, "basic", s2_ov, s1_ov, exp_ov,
                    csr_check="CHECK_CSRS_UNCHANGED_FP",
                )
                vd_init_fp = [_b32(100.0)] * 4
                mask_bits = 0b1010
                exp_masked = [
                    exp_ov[i] if (mask_bits >> i) & 1 else vd_init_fp[i]
                    for i in range(4)
                ]
                emit_binop_vv_masked(
                    tf, mnemonic, sew, "basic",
                    s2_ov, s1_ov, vd_init_fp, mask_bits, exp_masked,
                    csr_check="CHECK_CSRS_UNCHANGED_FP",
                )
        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP Binary VF ---
    for mnemonic, cfn, subdir in _FP_BINOP_VF:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / subdir / fname
        tf = TestFile(mnemonic, f"FP {mnemonic}")
        for sew in FP_SEWS:
            for name, s2, sc_bits in fp_binop_vf(sew):
                exp = [cfn(a, sc_bits, sew) for a in s2]
                emit_fp_binop_vf(tf, mnemonic, sew, name, s2, sc_bits, exp)

            # VF Overlap + Masked at e32 only
            if sew == 32:
                s2_ov = [_b32(1.0), _b32(2.0), _b32(3.0), _b32(4.0)]
                sc_ov = _b32(10.0)
                exp_ov = [cfn(a, sc_ov, sew) for a in s2_ov]
                emit_fp_binop_vf_overlap(
                    tf, mnemonic, sew, "basic", s2_ov, sc_ov, exp_ov,
                )
                vd_init_fp = [_b32(100.0)] * 4
                mask_bits = 0b1010
                exp_masked = [
                    exp_ov[i] if (mask_bits >> i) & 1 else vd_init_fp[i]
                    for i in range(4)
                ]
                emit_fp_binop_vf_masked(
                    tf, mnemonic, sew, "basic",
                    s2_ov, sc_ov, vd_init_fp, mask_bits, exp_masked,
                )
        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP FMA VV ---
    for mnemonic, cfn in _FP_FMA_VV:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_muladd" / fname
        tf = TestFile(mnemonic, f"FP FMA {mnemonic}")
        for sew in FP_SEWS:
            for name, vd_init, vs1, vs2 in fp_fma_vv(sew):
                exp = [cfn(d, s1, s2, sew) for d, s1, s2 in zip(vd_init, vs1, vs2)]
                emit_ternary_vv(tf, mnemonic, sew, name, vd_init, vs1, vs2, exp,
                                csr_check="CHECK_CSRS_UNCHANGED_FP")

            # Overlap (vd==vs2) + Masked at e32 only
            if sew == 32:
                vd_ov = [_b32(10.0), _b32(20.0), _b32(30.0), _b32(40.0)]
                vs1_ov = [_b32(1.0), _b32(2.0), _b32(3.0), _b32(4.0)]
                exp_ov = [cfn(d, s1, d, sew) for d, s1 in zip(vd_ov, vs1_ov)]
                emit_ternary_vv_overlap(
                    tf, mnemonic, sew, "basic", vd_ov, vs1_ov, exp_ov,
                    csr_check="CHECK_CSRS_UNCHANGED_FP",
                )
                vd_init_m = [_b32(100.0)] * 4
                vs2_m = [_b32(10.0), _b32(20.0), _b32(30.0), _b32(40.0)]
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
                    csr_check="CHECK_CSRS_UNCHANGED_FP",
                )
        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP Compare VV ---
    for mnemonic, cfn in _FP_CMP_VV:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_cmp" / fname
        tf = TestFile(mnemonic, f"FP compare {mnemonic}")
        for sew in FP_SEWS:
            for name, s2, s1 in fp_binop_vv(sew):
                mask = 0
                for i, (a, b) in enumerate(zip(s2, s1)):
                    if cfn(a, b, sew):
                        mask |= 1 << i
                emit_compare_vv(tf, mnemonic, sew, name, s2, s1, mask,
                                csr_check="CHECK_CSRS_UNCHANGED_FP")

            # Masked compare at e32 only
            if sew == 32:
                s2_m = [_b32(1.0), _b32(2.0), _b32(1.0), _b32(4.0)]
                s1_m = [_b32(1.0), _b32(1.0), _b32(2.0), _b32(4.0)]
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
                    csr_check="CHECK_CSRS_UNCHANGED_FP",
                )
        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP Compare VF ---
    for mnemonic, cfn in _FP_CMP_VF:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_cmp" / fname
        tf = TestFile(mnemonic, f"FP compare {mnemonic}")
        for sew in FP_SEWS:
            for name, s2, sc_bits in fp_binop_vf(sew):
                mask = 0
                for i, a in enumerate(s2):
                    if cfn(a, sc_bits, sew):
                        mask |= 1 << i
                emit_compare_vf(tf, mnemonic, sew, name, s2, sc_bits, mask)
        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP FMA VF ---
    for mnemonic, cfn in _FP_FMA_VF:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_muladd" / fname
        tf = TestFile(mnemonic, f"FP FMA {mnemonic}")
        for sew in FP_SEWS:
            for name, vd_init, sc_bits, vs2 in fp_fma_vf(sew):
                exp = [cfn(d, sc_bits, s2, sew) for d, s2 in zip(vd_init, vs2)]
                emit_ternary_vf(tf, mnemonic, sew, name, vd_init, sc_bits, vs2, exp,
                                csr_check="CHECK_CSRS_UNCHANGED_FP")
        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP sqrt ---
    from ..common import f32_to_bits, f64_to_bits
    for sew in FP_SEWS:
        mnemonic = "vfsqrt.v"
        fname = f"vfsqrt_e{sew}.S"
        fpath = base_dir / "tests" / "float_misc" / fname
        tf = TestFile(mnemonic, f"FP sqrt SEW={sew}")
        if sew == 32:
            src = [f32_to_bits(1.0), f32_to_bits(4.0), f32_to_bits(9.0), f32_to_bits(16.0)]
        else:
            src = [f64_to_bits(1.0), f64_to_bits(4.0), f64_to_bits(9.0), f64_to_bits(16.0)]
        exp = [fsqrt(v, sew) for v in src]
        nbytes = NUM_ELEMS * (sew // 8)
        cn_res = tf.next_check(f"vfsqrt.v e{sew}: result")
        tag = f"tc{cn_res}"

        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_src")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vfsqrt.v {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn_res}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED_FP")

        tf.data_align(sew)
        tf.data_label(f"{tag}_src", format_data_line(src, sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP classify ---
    for sew in FP_SEWS:
        import math
        mnemonic = "vfclass.v"
        fname = f"vfclass_e{sew}.S"
        fpath = base_dir / "tests" / "float_misc" / fname
        tf = TestFile(mnemonic, f"FP classify SEW={sew}")
        if sew == 32:
            b = f32_to_bits
            src = [b(1.0), b(-1.0), b(0.0), b(float('inf'))]
        else:
            b = f64_to_bits
            src = [b(1.0), b(-1.0), b(0.0), b(float('inf'))]
        exp = [fclass(v, sew) for v in src]
        nbytes = NUM_ELEMS * (sew // 8)
        cn_res = tf.next_check(f"vfclass.v e{sew}: result")
        tag = f"tc{cn_res}"

        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_src")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vfclass.v {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn_res}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED_FP")

        tf.data_align(sew)
        tf.data_label(f"{tag}_src", format_data_line(src, sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP reciprocal square-root estimate: vfrsqrt7.v ---
    # Use inputs whose 1/sqrt(x) is exactly representable as a float,
    # so the 7-bit-accurate hardware result matches our exact computation.
    # Powers of 4 give exact results: 1/sqrt(4^n) = 2^{-n}.
    from ..compute.floating_point import frsqrt7, frec7
    for sew in FP_SEWS:
        mnemonic = "vfrsqrt7.v"
        fname = f"vfrsqrt7_e{sew}.S"
        fpath = base_dir / "tests" / "float_misc" / fname
        tf = TestFile(mnemonic, f"FP reciprocal sqrt estimate SEW={sew}")
        if sew == 32:
            # 1/sqrt(1)=1, 1/sqrt(4)=0.5, 1/sqrt(0.25)=2, 1/sqrt(16)=0.25
            src = [f32_to_bits(1.0), f32_to_bits(4.0),
                   f32_to_bits(0.25), f32_to_bits(16.0)]
        else:
            src = [f64_to_bits(1.0), f64_to_bits(4.0),
                   f64_to_bits(0.25), f64_to_bits(16.0)]
        exp = [frsqrt7(v, sew) for v in src]
        nbytes = NUM_ELEMS * (sew // 8)
        cn_res = tf.next_check(f"vfrsqrt7.v e{sew}: result")
        tag = f"tc{cn_res}"

        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_src")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vfrsqrt7.v {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn_res}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED_FP")

        tf.data_align(sew)
        tf.data_label(f"{tag}_src", format_data_line(src, sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP reciprocal estimate: vfrec7.v ---
    # Use powers of 2 so 1/x is exact: 1/1=1, 1/2=0.5, 1/4=0.25, 1/0.5=2.
    for sew in FP_SEWS:
        mnemonic = "vfrec7.v"
        fname = f"vfrec7_e{sew}.S"
        fpath = base_dir / "tests" / "float_misc" / fname
        tf = TestFile(mnemonic, f"FP reciprocal estimate SEW={sew}")
        if sew == 32:
            src = [f32_to_bits(1.0), f32_to_bits(2.0),
                   f32_to_bits(4.0), f32_to_bits(0.5)]
        else:
            src = [f64_to_bits(1.0), f64_to_bits(2.0),
                   f64_to_bits(4.0), f64_to_bits(0.5)]
        exp = [frec7(v, sew) for v in src]
        nbytes = NUM_ELEMS * (sew // 8)
        cn_res = tf.next_check(f"vfrec7.v e{sew}: result")
        tag = f"tc{cn_res}"

        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_src")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vfrec7.v {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn_res}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED_FP")

        tf.data_align(sew)
        tf.data_label(f"{tag}_src", format_data_line(src, sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP conversions (single-width) ---
    _CVT_OPS: list[tuple[str, Callable]] = [
        ("vfcvt.xu.f.v",     fcvt_xu_f),
        ("vfcvt.x.f.v",      fcvt_x_f),
        ("vfcvt.rtz.xu.f.v", fcvt_rtz_xu_f),
        ("vfcvt.rtz.x.f.v",  fcvt_rtz_x_f),
        ("vfcvt.f.xu.v",     fcvt_f_xu),
        ("vfcvt.f.x.v",      fcvt_f_x),
    ]
    for mnemonic, cfn in _CVT_OPS:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_convert" / fname
        tf = TestFile(mnemonic, f"FP conversion {mnemonic}")
        for sew in FP_SEWS:
            if "f.xu" in mnemonic or "f.x" in mnemonic:
                # Input is integer
                src = [0, 1, 42, (1 << (sew - 1)) - 1]
            else:
                # Input is FP
                if sew == 32:
                    b = f32_to_bits
                    src = [b(0.0), b(1.0), b(42.5), b(-3.75)]
                else:
                    b = f64_to_bits
                    src = [b(0.0), b(1.0), b(42.5), b(-3.75)]
            exp = [cfn(v, sew) for v in src]
            nbytes = NUM_ELEMS * (sew // 8)
            cn_res = tf.next_check(f"{mnemonic} e{sew}: result")
            tag = f"tc{cn_res}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}")
            tf.code(f"SET_TEST_NUM {cn_res}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED_FP")

            tf.data_align(sew)
            tf.data_label(f"{tag}_src", format_data_line(src, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        tf.write(fpath)
        generated.append(str(fpath))

    return generated
