from __future__ import annotations

"""Generate tests for integer arithmetic instructions:
vadd, vsub, vrsub, vand, vor, vxor, vsll, vsrl, vsra,
vminu, vmin, vmaxu, vmax, vmul, vmulh, vmulhu, vmulhsu,
vdivu, vdiv, vremu, vrem.
"""

from pathlib import Path
from typing import Callable

from ..common import SEWS, M, U
from ..testfile import TestFile
from ..emit import (
    emit_binop_vv, emit_binop_vx, emit_binop_vi,
    emit_binop_vv_overlap, emit_binop_vv_masked,
    emit_binop_vv_overlap_vs1, emit_binop_vv_self_overlap,
    emit_binop_vx_overlap, emit_binop_vx_masked,
    emit_binop_vi_overlap, emit_binop_vi_masked,
)
from ..vectors import binop_vv, binop_vx, binop_vi
from ..compute.integer import (
    add, sub, rsub, and_, or_, xor, sll, srl, sra,
    minu, min_, maxu, max_, mul, mulh, mulhu, mulhsu,
    divu, div, remu, rem,
)

# (base_mnemonic, forms, compute_fn, subdirectory, flags)
_INSTRUCTIONS: list[tuple[str, list[str], Callable, str, dict]] = [
    # Arithmetic
    ("vadd",    ["vv", "vx", "vi"], add,    "int_arith", {}),
    ("vsub",    ["vv", "vx"],       sub,    "int_arith", {}),
    ("vrsub",   ["vx", "vi"],       rsub,   "int_arith", {}),
    # Logical
    ("vand",    ["vv", "vx", "vi"], and_,   "int_logical", {}),
    ("vor",     ["vv", "vx", "vi"], or_,    "int_logical", {}),
    ("vxor",    ["vv", "vx", "vi"], xor,    "int_logical", {}),
    # Shift
    ("vsll",    ["vv", "vx", "vi"], sll,    "int_shift", {"unsigned_imm": True}),
    ("vsrl",    ["vv", "vx", "vi"], srl,    "int_shift", {"unsigned_imm": True}),
    ("vsra",    ["vv", "vx", "vi"], sra,    "int_shift", {"unsigned_imm": True}),
    # Min / Max
    ("vminu",   ["vv", "vx"],       minu,   "int_minmax", {}),
    ("vmin",    ["vv", "vx"],       min_,   "int_minmax", {}),
    ("vmaxu",   ["vv", "vx"],       maxu,   "int_minmax", {}),
    ("vmax",    ["vv", "vx"],       max_,   "int_minmax", {}),
    # Multiply
    ("vmul",    ["vv", "vx"],       mul,    "int_mul", {}),
    ("vmulh",   ["vv", "vx"],       mulh,   "int_mul", {}),
    ("vmulhu",  ["vv", "vx"],       mulhu,  "int_mul", {}),
    ("vmulhsu", ["vv", "vx"],       mulhsu, "int_mul", {}),
    # Divide
    ("vdivu",   ["vv", "vx"],       divu,   "int_div", {}),
    ("vdiv",    ["vv", "vx"],       div,    "int_div", {}),
    ("vremu",   ["vv", "vx"],       remu,   "int_div", {}),
    ("vrem",    ["vv", "vx"],       rem,    "int_div", {}),
]


def generate(base_dir: Path) -> list[str]:
    """Generate all integer arithmetic test files. Returns list of paths."""
    generated: list[str] = []

    for base, forms, cfn, subdir, flags in _INSTRUCTIONS:
        for form in forms:
            mnemonic = f"{base}.{form}"
            fname = f"{base}_{form}.S"
            fpath = base_dir / "tests" / subdir / fname

            tf = TestFile(mnemonic, f"Integer {base} ({form} form)")

            for sew in SEWS:
                if form == "vv":
                    for name, s2, s1 in binop_vv(sew):
                        exp = [cfn(a, b, sew) for a, b in zip(s2, s1)]
                        emit_binop_vv(tf, mnemonic, sew, name, s2, s1, exp)

                    # Overlap test (vd==vs2) at e32 only
                    if sew == 32:
                        s2_ov, s1_ov = [1, 2, 3, 4], [10, 20, 30, 40]
                        exp_ov = [cfn(a, b, sew) for a, b in zip(s2_ov, s1_ov)]
                        emit_binop_vv_overlap(
                            tf, mnemonic, sew, "basic", s2_ov, s1_ov, exp_ov,
                        )

                        # Masked test (v0.t, mask=0b1010) at e32 only
                        vd_init = [0xDEAD] * 4
                        mask_bits = 0b1010  # elements 1,3 active
                        exp_m = [cfn(a, b, sew) for a, b in zip(s2_ov, s1_ov)]
                        exp_masked = [
                            exp_m[i] if (mask_bits >> i) & 1 else vd_init[i]
                            for i in range(4)
                        ]
                        emit_binop_vv_masked(
                            tf, mnemonic, sew, "basic",
                            s2_ov, s1_ov, vd_init, mask_bits, exp_masked,
                        )

                        # Overlap test (vd==vs1) at e32
                        exp_vs1 = [cfn(a, b, sew) for a, b in zip(s2_ov, s1_ov)]
                        emit_binop_vv_overlap_vs1(
                            tf, mnemonic, sew, "basic", s2_ov, s1_ov, exp_vs1,
                        )

                        # Self-overlap (vd==vs2==vs1) at e32
                        src_self = [1, 2, 3, 4]
                        exp_self = [cfn(a, a, sew) for a in src_self]
                        emit_binop_vv_self_overlap(
                            tf, mnemonic, sew, "basic", src_self, exp_self,
                        )

                elif form == "vx":
                    for name, vec, sc in binop_vx(sew):
                        exp = [cfn(a, sc, sew) for a in vec]
                        emit_binop_vx(tf, mnemonic, sew, name, vec, U(sc, sew), exp)

                    # Overlap + Masked at e32 only
                    if sew == 32:
                        s2_ov = [1, 2, 3, 4]
                        sc_ov = U(10, sew)
                        exp_ov = [cfn(a, sc_ov, sew) for a in s2_ov]
                        emit_binop_vx_overlap(
                            tf, mnemonic, sew, "basic", s2_ov, sc_ov, exp_ov,
                        )
                        vd_init = [0xDEAD] * 4
                        mask_bits = 0b1010
                        exp_masked = [
                            exp_ov[i] if (mask_bits >> i) & 1 else vd_init[i]
                            for i in range(4)
                        ]
                        emit_binop_vx_masked(
                            tf, mnemonic, sew, "basic",
                            s2_ov, sc_ov, vd_init, mask_bits, exp_masked,
                        )

                elif form == "vi":
                    uimm = flags.get("unsigned_imm", False)
                    for name, vec, imm in binop_vi(sew, unsigned_imm=uimm):
                        imm_val = U(imm, sew)
                        exp = [cfn(a, imm_val, sew) for a in vec]
                        emit_binop_vi(tf, mnemonic, sew, name, vec, imm, exp)

                    # Overlap + Masked at e32 only
                    if sew == 32:
                        s2_ov = [1, 2, 3, 4]
                        imm_ov = 1
                        imm_val_ov = U(imm_ov, sew)
                        exp_ov = [cfn(a, imm_val_ov, sew) for a in s2_ov]
                        emit_binop_vi_overlap(
                            tf, mnemonic, sew, "basic", s2_ov, imm_ov, exp_ov,
                        )
                        vd_init = [0xDEAD] * 4
                        mask_bits = 0b1010
                        exp_masked = [
                            exp_ov[i] if (mask_bits >> i) & 1 else vd_init[i]
                            for i in range(4)
                        ]
                        emit_binop_vi_masked(
                            tf, mnemonic, sew, "basic",
                            s2_ov, imm_ov, vd_init, mask_bits, exp_masked,
                        )

            tf.write(fpath)
            generated.append(str(fpath))

    return generated
