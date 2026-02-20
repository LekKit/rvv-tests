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
from ..emit import emit_binop_vv, emit_binop_vx, emit_binop_vi
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

                elif form == "vx":
                    for name, vec, sc in binop_vx(sew):
                        exp = [cfn(a, sc, sew) for a in vec]
                        emit_binop_vx(tf, mnemonic, sew, name, vec, U(sc, sew), exp)

                elif form == "vi":
                    uimm = flags.get("unsigned_imm", False)
                    for name, vec, imm in binop_vi(sew, unsigned_imm=uimm):
                        # For signed imm, the HW sign-extends the 5-bit value
                        imm_val = U(imm, sew) if imm >= 0 else U(imm, sew)
                        exp = [cfn(a, imm_val, sew) for a in vec]
                        emit_binop_vi(tf, mnemonic, sew, name, vec, imm, exp)

            tf.write(fpath)
            generated.append(str(fpath))

    return generated
