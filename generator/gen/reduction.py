from __future__ import annotations

"""Generate tests for reduction instructions:
vredsum, vredmaxu, vredmax, vredminu, vredmin, vredand, vredor, vredxor,
vwredsumu, vwredsum,
vfredosum, vfredusum, vfredmax, vfredmin,
vfwredosum, vfwredusum.
"""

from pathlib import Path
from functools import reduce
from typing import Callable

from ..common import SEWS, FP_SEWS, WIDENING_SEWS, M, S, U, NUM_ELEMS, format_data_line
from ..testfile import TestFile
from ..emit import emit_reduction, VREG_DST, VREG_SRC2, VREG_SRC1
from ..vectors import reduction_int
from ..compute.floating_point import fredosum, fredusum, fredmax, fredmin


def _ireduce(op: Callable, init: int, vec: list[int], sew: int) -> int:
    """Compute integer reduction."""
    acc = init
    for v in vec:
        acc = op(acc, v, sew)
    return U(acc, sew)


def _isum(a: int, b: int, sew: int) -> int:
    return U(a + b, sew)

def _imaxu(a: int, b: int, sew: int) -> int:
    return max(U(a, sew), U(b, sew))

def _imax(a: int, b: int, sew: int) -> int:
    return U(max(S(a, sew), S(b, sew)), sew)

def _iminu(a: int, b: int, sew: int) -> int:
    return min(U(a, sew), U(b, sew))

def _imin(a: int, b: int, sew: int) -> int:
    return U(min(S(a, sew), S(b, sew)), sew)

def _iand(a: int, b: int, sew: int) -> int:
    return a & b & M(sew)

def _ior(a: int, b: int, sew: int) -> int:
    return (a | b) & M(sew)

def _ixor(a: int, b: int, sew: int) -> int:
    return (a ^ b) & M(sew)


_INT_REDUCTIONS: list[tuple[str, Callable]] = [
    ("vredsum.vs",  _isum),
    ("vredmaxu.vs", _imaxu),
    ("vredmax.vs",  _imax),
    ("vredminu.vs", _iminu),
    ("vredmin.vs",  _imin),
    ("vredand.vs",  _iand),
    ("vredor.vs",   _ior),
    ("vredxor.vs",  _ixor),
]


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "reduction"

    # --- Integer reductions ---
    for mnemonic, op_fn in _INT_REDUCTIONS:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Reduction {mnemonic}")
        for sew in SEWS:
            for name, init, vec in reduction_int(sew):
                exp = _ireduce(op_fn, init, vec, sew)
                emit_reduction(tf, mnemonic, sew, name, init, vec, exp)
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- FP reductions ---
    from ..common import f32_to_bits, f64_to_bits
    _FP_RED: list[tuple[str, Callable]] = [
        ("vfredosum.vs", fredosum),
        ("vfredusum.vs", fredusum),
        ("vfredmax.vs",  fredmax),
        ("vfredmin.vs",  fredmin),
    ]
    for mnemonic, cfn in _FP_RED:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"FP reduction {mnemonic}")
        for sew in FP_SEWS:
            if sew == 32:
                b = f32_to_bits
                init = b(0.0)
                vec = [b(1.0), b(2.0), b(3.0), b(4.0)]
            else:
                b = f64_to_bits
                init = b(0.0)
                vec = [b(1.0), b(2.0), b(3.0), b(4.0)]
            exp = cfn(init, vec, sew)
            emit_reduction(tf, mnemonic, sew, "basic", init, vec, exp,
                           csr_check="CHECK_CSRS_UNCHANGED_FP")
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Widening integer reductions ---
    # vwredsumu.vs: unsigned widening sum (source SEW, accumulator 2*SEW)
    # vwredsum.vs:  signed widening sum
    for mnemonic, is_signed in [("vwredsumu.vs", False), ("vwredsum.vs", True)]:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Widening reduction {mnemonic}")
        for sew in WIDENING_SEWS:
            dsew = 2 * sew
            m = M(sew)
            # Test data: source elements at SEW, accumulator at 2*SEW
            test_cases_w = [
                ("basic", 0, [1, 2, 3, 4]),
                ("init",  100, [1, 2, 3, 4]),
                ("max",   0, [m, m, m, m]),
            ]
            for name, init_val, vec in test_cases_w:
                # Compute widening sum
                acc = U(init_val, dsew)
                for v in vec:
                    if is_signed:
                        acc = U(S(acc, dsew) + S(v, sew), dsew)
                    else:
                        acc = U(U(acc, dsew) + U(v, sew), dsew)
                exp_scalar = acc

                nbytes_d = dsew // 8  # check element 0 only
                cn = tf.next_check(f"{mnemonic} e{sew} {name}: result")
                tag = f"tc{cn}"

                tf.blank()
                tf.comment(f"Test {cn}: {mnemonic} SEW={sew} {name}")
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                # Load source vector
                tf.code(f"la t1, {tag}_vec")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                # Load init scalar into vs1 at dest width
                tf.code(f"vsetvli t0, t0, e{dsew}, m1, tu, mu")
                tf.code(f"la t1, {tag}_init")
                tf.code(f"vle{dsew}.v {VREG_SRC1}, (t1)")
                # Restore source vtype for the instruction
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
                tf.code(f"SET_TEST_NUM {cn}")
                # Read result at dest width
                tf.code(f"vsetvli t0, t0, e{dsew}, m1, tu, mu")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes_d}")
                tf.code("CHECK_VSTART_ZERO")

                tf.data_align(dsew)
                tf.data_label(f"{tag}_vec", format_data_line(vec, sew))
                tf.data_label(f"{tag}_init",
                              format_data_line(
                                  [init_val] + [0] * (NUM_ELEMS - 1), dsew))
                tf.data_label(f"{tag}_exp",
                              format_data_line(
                                  [exp_scalar] + [0] * (NUM_ELEMS - 1), dsew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- Widening FP reductions ---
    from ..compute.floating_point import fwredosum, fwredusum
    _WFP_RED: list[tuple[str, object]] = [
        ("vfwredosum.vs", fwredosum),
        ("vfwredusum.vs", fwredusum),
    ]
    for mnemonic, cfn in _WFP_RED:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Widening FP reduction {mnemonic}")
        sew = 32
        dsew = 64
        b32 = f32_to_bits
        b64 = f64_to_bits
        vec = [b32(1.0), b32(2.0), b32(3.0), b32(4.0)]
        init_val = b64(0.0)
        exp = cfn(init_val, vec, sew)

        cn = tf.next_check(f"{mnemonic}: result")
        tag = f"tc{cn}"

        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_vec")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        # Load init at double width
        tf.code(f"vsetvli t0, t0, e{dsew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_init")
        tf.code(f"vle{dsew}.v {VREG_SRC1}, (t1)")
        # Restore source vtype
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code("SAVE_CSRS")
        tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"vsetvli t0, t0, e{dsew}, m1, tu, mu")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {dsew // 8}")
        tf.code("CHECK_VSTART_ZERO")
        tf.code("CHECK_FCSR_UNCHANGED")

        tf.data_align(dsew)
        tf.data_label(f"{tag}_vec", format_data_line(vec, sew))
        tf.data_label(f"{tag}_init",
                      format_data_line([init_val] + [0] * (NUM_ELEMS - 1), dsew))
        tf.data_label(f"{tag}_exp",
                      format_data_line([exp] + [0] * (NUM_ELEMS - 1), dsew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    return generated
