from __future__ import annotations

"""Generate tests for widening FP instructions:
vfwadd.vv/.vf/.wv/.wf, vfwsub.vv/.vf/.wv/.wf, vfwmul.vv/.vf,
vfwmacc.vv/.vf, vfwnmacc.vv/.vf, vfwmsac.vv/.vf, vfwnmsac.vv/.vf,
vfwcvt.*, vfncvt.*.
"""

from pathlib import Path
from typing import Callable

from ..common import (
    NUM_ELEMS, M, S, U,
    format_data_line, format_hex,
    f32_to_bits, f64_to_bits, bits_to_f32, bits_to_f64,
)
from ..testfile import TestFile
from ..emit import VREG_DST, VREG_SRC2, VREG_SRC1
from ..compute.floating_point import (
    fwadd, fwsub, fwmul, fwadd_w, fwsub_w,
    fwmacc, fwnmacc, fwmsac, fwnmsac,
    fwcvt_f_f, fwcvt_xu_f, fwcvt_x_f,
    fwcvt_rtz_xu_f, fwcvt_rtz_x_f, fwcvt_f_xu, fwcvt_f_x,
    fncvt_f_f, fncvt_xu_f, fncvt_x_f,
    fncvt_rtz_xu_f, fncvt_rtz_x_f, fncvt_f_xu, fncvt_f_x,
    fncvt_rod_f_f,
)

# Only SEW=32 can widen to f64.
_SRC_SEW = 32
_DST_SEW = 64


def _fp32_test_vectors() -> list[tuple[str, list[int], list[int]]]:
    b = f32_to_bits
    return [
        ("basic",  [b(1.0), b(2.0), b(3.0), b(4.0)],
                   [b(0.5), b(0.25), b(0.125), b(0.0625)]),
        ("negpos", [b(-1.0), b(-2.0), b(3.0), b(4.0)],
                   [b(1.0), b(2.0), b(-3.0), b(-4.0)]),
    ]


def _fp32_vf_vectors() -> list[tuple[str, list[int], int]]:
    b = f32_to_bits
    return [
        ("basic", [b(1.0), b(2.0), b(3.0), b(4.0)], b(10.0)),
        ("neg",   [b(1.0), b(2.0), b(3.0), b(4.0)], b(-1.0)),
    ]


def _fp32_fma_vectors() -> list[tuple[str, list[int], list[int], list[int]]]:
    """(name, vd_init_f64_bits[4], vs1_f32_bits[4], vs2_f32_bits[4])."""
    b32 = f32_to_bits
    b64 = f64_to_bits
    return [
        ("basic",
         [b64(0.0)] * 4,
         [b32(1.0), b32(2.0), b32(3.0), b32(4.0)],
         [b32(10.0), b32(20.0), b32(30.0), b32(40.0)]),
        ("accum",
         [b64(100.0), b64(200.0), b64(300.0), b64(400.0)],
         [b32(1.0), b32(1.0), b32(1.0), b32(1.0)],
         [b32(1.0), b32(2.0), b32(3.0), b32(4.0)]),
    ]


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []

    sew = _SRC_SEW
    dsew = _DST_SEW

    # ------------------------------------------------------------------
    # Widening FP arith VV
    # ------------------------------------------------------------------
    _WFOP_VV: list[tuple[str, Callable]] = [
        ("vfwadd.vv", fwadd),
        ("vfwsub.vv", fwsub),
        ("vfwmul.vv", fwmul),
    ]
    for mnemonic, cfn in _WFOP_VV:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_widening" / fname
        tf = TestFile(mnemonic, f"Widening FP {mnemonic}")
        for name, s2, s1 in _fp32_test_vectors():
            exp = [cfn(a, b, sew) for a, b in zip(s2, s1)]
            nbytes = NUM_ELEMS * (dsew // 8)
            cn = tf.next_check(f"{mnemonic} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_VSTART_ZERO")
            tf.code("CHECK_FCSR_UNCHANGED")

            tf.data_align(dsew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))
        tf.write(fpath)
        generated.append(str(fpath))

    # ------------------------------------------------------------------
    # Widening FP arith VF
    # ------------------------------------------------------------------
    _WFOP_VF: list[tuple[str, Callable]] = [
        ("vfwadd.vf", fwadd),
        ("vfwsub.vf", fwsub),
        ("vfwmul.vf", fwmul),
    ]
    for mnemonic, cfn in _WFOP_VF:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_widening" / fname
        tf = TestFile(mnemonic, f"Widening FP {mnemonic}")
        for name, s2, sc_bits in _fp32_vf_vectors():
            exp = [cfn(a, sc_bits, sew) for a in s2]
            nbytes = NUM_ELEMS * (dsew // 8)
            cn = tf.next_check(f"{mnemonic} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_fsc")
            tf.code(f"flw fa0, 0(t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, fa0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_VSTART_ZERO")
            tf.code("CHECK_FCSR_UNCHANGED")

            tf.data_align(dsew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_fsc", format_data_line([sc_bits], sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))
        tf.write(fpath)
        generated.append(str(fpath))

    # ------------------------------------------------------------------
    # Widening FP arith .WV (vs2 is 2*SEW, vs1 is SEW)
    # ------------------------------------------------------------------
    _WFOP_WV: list[tuple[str, Callable]] = [
        ("vfwadd.wv", fwadd_w),
        ("vfwsub.wv", fwsub_w),
    ]
    b32 = f32_to_bits
    b64 = f64_to_bits
    wv_data = [
        ("basic",
         [b64(10.0), b64(20.0), b64(30.0), b64(40.0)],
         [b32(1.0), b32(2.0), b32(3.0), b32(4.0)]),
        ("neg",
         [b64(-1.0), b64(-2.0), b64(3.0), b64(4.0)],
         [b32(1.0), b32(2.0), b32(-3.0), b32(-4.0)]),
    ]
    for mnemonic, cfn in _WFOP_WV:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_widening" / fname
        tf = TestFile(mnemonic, f"Widening FP {mnemonic}")
        for name, s2_wide, s1 in wv_data:
            exp = [cfn(a, b, sew) for a, b in zip(s2_wide, s1)]
            nbytes = NUM_ELEMS * (dsew // 8)
            cn = tf.next_check(f"{mnemonic} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            # Load wide source
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, {tag}_s2w")
            tf.code(f"vle{dsew}.v {VREG_SRC2}, (t1)")
            # Load narrow source
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_VSTART_ZERO")
            tf.code("CHECK_FCSR_UNCHANGED")

            tf.data_align(dsew)
            tf.data_label(f"{tag}_s2w", format_data_line(s2_wide, dsew))
            tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))
        tf.write(fpath)
        generated.append(str(fpath))

    # ------------------------------------------------------------------
    # Widening FP FMA VV
    # ------------------------------------------------------------------
    _WFMA_VV: list[tuple[str, Callable]] = [
        ("vfwmacc.vv",  fwmacc),
        ("vfwnmacc.vv", fwnmacc),
        ("vfwmsac.vv",  fwmsac),
        ("vfwnmsac.vv", fwnmsac),
    ]
    for mnemonic, cfn in _WFMA_VV:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_widening" / fname
        tf = TestFile(mnemonic, f"Widening FP FMA {mnemonic}")
        for name, vd_init, vs1, vs2 in _fp32_fma_vectors():
            exp = [cfn(d, s1, s2, sew) for d, s1, s2 in zip(vd_init, vs1, vs2)]
            nbytes = NUM_ELEMS * (dsew // 8)
            cn = tf.next_check(f"{mnemonic} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            # Load vd (wide)
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, {tag}_vd")
            tf.code(f"vle{dsew}.v {VREG_DST}, (t1)")
            # Load vs1, vs2 (narrow)
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC1}, {VREG_SRC2}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_VSTART_ZERO")
            tf.code("CHECK_FCSR_UNCHANGED")

            tf.data_align(dsew)
            tf.data_label(f"{tag}_vd", format_data_line(vd_init, dsew))
            tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
            tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))
        tf.write(fpath)
        generated.append(str(fpath))

    # ------------------------------------------------------------------
    # Widening FP FMA VF
    # ------------------------------------------------------------------
    _WFMA_VF: list[tuple[str, Callable]] = [
        ("vfwmacc.vf",  fwmacc),
        ("vfwnmacc.vf", fwnmacc),
        ("vfwmsac.vf",  fwmsac),
        ("vfwnmsac.vf", fwnmsac),
    ]
    _fp32_fma_vf_data: list[tuple[str, list[int], int, list[int]]] = [
        ("basic",
         [b64(0.0)] * 4,
         b32(2.0),
         [b32(10.0), b32(20.0), b32(30.0), b32(40.0)]),
        ("accum",
         [b64(100.0), b64(200.0), b64(300.0), b64(400.0)],
         b32(1.0),
         [b32(1.0), b32(2.0), b32(3.0), b32(4.0)]),
    ]
    for mnemonic, cfn in _WFMA_VF:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_widening" / fname
        tf = TestFile(mnemonic, f"Widening FP FMA {mnemonic}")
        for name, vd_init, sc_bits, vs2 in _fp32_fma_vf_data:
            exp = [cfn(d, sc_bits, s2, sew) for d, s2 in zip(vd_init, vs2)]
            nbytes = NUM_ELEMS * (dsew // 8)
            cn = tf.next_check(f"{mnemonic} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            # Load vd (wide)
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, {tag}_vd")
            tf.code(f"vle{dsew}.v {VREG_DST}, (t1)")
            # Load vs2 (narrow)
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            # Load FP scalar
            tf.code(f"la t1, {tag}_fsc")
            tf.code(f"flw fa0, 0(t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, fa0, {VREG_SRC2}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_VSTART_ZERO")
            tf.code("CHECK_FCSR_UNCHANGED")

            tf.data_align(dsew)
            tf.data_label(f"{tag}_vd", format_data_line(vd_init, dsew))
            tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
            tf.data_label(f"{tag}_fsc", format_data_line([sc_bits], sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))
        tf.write(fpath)
        generated.append(str(fpath))

    # ------------------------------------------------------------------
    # Widening FP arith .WF (vs2 is 2*SEW, scalar is SEW)
    # ------------------------------------------------------------------
    _WFOP_WF: list[tuple[str, Callable]] = [
        ("vfwadd.wf", fwadd_w),
        ("vfwsub.wf", fwsub_w),
    ]
    wf_data: list[tuple[str, list[int], int]] = [
        ("basic",
         [b64(10.0), b64(20.0), b64(30.0), b64(40.0)],
         b32(1.0)),
        ("neg",
         [b64(-1.0), b64(-2.0), b64(3.0), b64(4.0)],
         b32(-3.0)),
    ]
    for mnemonic, cfn in _WFOP_WF:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_widening" / fname
        tf = TestFile(mnemonic, f"Widening FP {mnemonic}")
        for name, s2_wide, sc_bits in wf_data:
            exp = [cfn(a, sc_bits, sew) for a in s2_wide]
            nbytes = NUM_ELEMS * (dsew // 8)
            cn = tf.next_check(f"{mnemonic} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            # Load wide source
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, {tag}_s2w")
            tf.code(f"vle{dsew}.v {VREG_SRC2}, (t1)")
            # Set vtype to narrow SEW for the instruction
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            # Load FP scalar
            tf.code(f"la t1, {tag}_fsc")
            tf.code(f"flw fa0, 0(t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, fa0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_VSTART_ZERO")
            tf.code("CHECK_FCSR_UNCHANGED")

            tf.data_align(dsew)
            tf.data_label(f"{tag}_s2w", format_data_line(s2_wide, dsew))
            tf.data_label(f"{tag}_fsc", format_data_line([sc_bits], sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))
        tf.write(fpath)
        generated.append(str(fpath))

    # ------------------------------------------------------------------
    # Widening FP conversions  (source f32 → dest f64 or int64)
    # ------------------------------------------------------------------
    _WFCVT: list[tuple[str, Callable, bool]] = [
        # (mnemonic, compute_fn, source_is_fp)
        ("vfwcvt.f.f.v",       fwcvt_f_f,       True),
        ("vfwcvt.xu.f.v",      fwcvt_xu_f,      True),
        ("vfwcvt.x.f.v",       fwcvt_x_f,       True),
        ("vfwcvt.rtz.xu.f.v",  fwcvt_rtz_xu_f,  True),
        ("vfwcvt.rtz.x.f.v",   fwcvt_rtz_x_f,   True),
        ("vfwcvt.f.xu.v",      fwcvt_f_xu,      False),
        ("vfwcvt.f.x.v",       fwcvt_f_x,       False),
    ]
    for mnemonic, cfn, src_is_fp in _WFCVT:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_widening" / fname
        tf = TestFile(mnemonic, f"Widening FP conversion {mnemonic}")

        if src_is_fp:
            # Use exact integer-valued floats to avoid INEXACT/NV flags
            src = [b32(0.0), b32(1.0), b32(42.0), b32(100.0)]
        else:
            src = [0, 1, 42, (1 << (sew - 1)) - 1]
        exp = [cfn(v, sew) for v in src]
        nbytes = NUM_ELEMS * (dsew // 8)
        cn = tf.next_check(f"{mnemonic}: result")
        tag = f"tc{cn}"

        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_src")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_VSTART_ZERO")
        tf.code("CHECK_FCSR_UNCHANGED")

        tf.data_align(dsew)
        tf.data_label(f"{tag}_src", format_data_line(src, sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))
        tf.write(fpath)
        generated.append(str(fpath))

    # ------------------------------------------------------------------
    # Narrowing FP conversions  (source f64/int64 → dest f32/int32)
    # ------------------------------------------------------------------
    _NFCVT: list[tuple[str, Callable, bool]] = [
        ("vfncvt.f.f.w",       fncvt_f_f,       True),
        ("vfncvt.xu.f.w",      fncvt_xu_f,      True),
        ("vfncvt.x.f.w",       fncvt_x_f,       True),
        ("vfncvt.rtz.xu.f.w",  fncvt_rtz_xu_f,  True),
        ("vfncvt.rtz.x.f.w",   fncvt_rtz_x_f,   True),
        ("vfncvt.f.xu.w",      fncvt_f_xu,      False),
        ("vfncvt.f.x.w",       fncvt_f_x,       False),
        ("vfncvt.rod.f.f.w",   fncvt_rod_f_f,   True),
    ]
    for mnemonic, cfn, src_is_fp in _NFCVT:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = base_dir / "tests" / "float_narrowing" / fname
        tf = TestFile(mnemonic, f"Narrowing FP conversion {mnemonic}")

        if src_is_fp:
            # Use exact integer-valued floats to avoid INEXACT/NV flags
            src = [b64(0.0), b64(1.0), b64(42.0), b64(100.0)]
        else:
            # Use values exactly representable in f32 to avoid INEXACT
            src = [0, 1, 42, 1000]
        exp = [cfn(v, sew) for v in src]
        # dest is 32-bit
        nbytes = NUM_ELEMS * (sew // 8)
        cn = tf.next_check(f"{mnemonic}: result")
        tag = f"tc{cn}"

        # Load source at 2*SEW = 64
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{dsew}, m2, tu, mu")
        tf.code(f"la t1, {tag}_src")
        tf.code(f"vle{dsew}.v {VREG_SRC2}, (t1)")
        # Set vtype to dest SEW for the narrowing instruction
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code("SAVE_CSRS")
        tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_VSTART_ZERO")
        tf.code("CHECK_FCSR_UNCHANGED")

        tf.data_align(dsew)
        tf.data_label(f"{tag}_src", format_data_line(src, dsew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
        tf.write(fpath)
        generated.append(str(fpath))

    return generated
