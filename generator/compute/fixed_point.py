from __future__ import annotations

"""Compute functions for RVV fixed-point arithmetic instructions."""

from ..common import M, S, U, roundoff_unsigned, roundoff_signed

# ===================================================================
# Saturating add / subtract  →  (result, vxsat_bit)
# ===================================================================

def saddu(a: int, b: int, sew: int) -> tuple[int, int]:
    r = U(a, sew) + U(b, sew)
    if r > M(sew):
        return M(sew), 1
    return r, 0

def sadd(a: int, b: int, sew: int) -> tuple[int, int]:
    r = S(a, sew) + S(b, sew)
    smax = (1 << (sew - 1)) - 1
    smin = -(1 << (sew - 1))
    if r > smax:
        return U(smax, sew), 1
    if r < smin:
        return U(smin, sew), 1
    return U(r, sew), 0

def ssubu(a: int, b: int, sew: int) -> tuple[int, int]:
    r = U(a, sew) - U(b, sew)
    if r < 0:
        return 0, 1
    return U(r, sew), 0

def ssub(a: int, b: int, sew: int) -> tuple[int, int]:
    r = S(a, sew) - S(b, sew)
    smax = (1 << (sew - 1)) - 1
    smin = -(1 << (sew - 1))
    if r > smax:
        return U(smax, sew), 1
    if r < smin:
        return U(smin, sew), 1
    return U(r, sew), 0

# ===================================================================
# Averaging add / subtract  (uses vxrm rounding mode)
# ===================================================================

def aaddu(a: int, b: int, sew: int, vxrm: int = 0) -> int:
    return U(roundoff_unsigned(U(a, sew) + U(b, sew), 1, vxrm), sew)

def aadd(a: int, b: int, sew: int, vxrm: int = 0) -> int:
    return U(roundoff_signed(S(a, sew) + S(b, sew), 1, vxrm), sew)

def asubu(a: int, b: int, sew: int, vxrm: int = 0) -> int:
    return U(roundoff_signed(U(a, sew) - U(b, sew), 1, vxrm), sew)

def asub(a: int, b: int, sew: int, vxrm: int = 0) -> int:
    return U(roundoff_signed(S(a, sew) - S(b, sew), 1, vxrm), sew)

# ===================================================================
# Fractional multiply with rounding and saturation  →  (result, vxsat)
# ===================================================================

def smul(a: int, b: int, sew: int, vxrm: int = 0) -> tuple[int, int]:
    """vsmul: signed fractional multiply.

    result = roundoff_signed(a * b, SEW-1, vxrm), then saturate.
    Special case: if both inputs are the minimum signed value, result
    saturates to the maximum positive value and vxsat is set.
    """
    sa, sb = S(a, sew), S(b, sew)
    smin = -(1 << (sew - 1))
    smax = (1 << (sew - 1)) - 1
    if sa == smin and sb == smin:
        return U(smax, sew), 1
    product = sa * sb
    r = roundoff_signed(product, sew - 1, vxrm)
    if r > smax:
        return U(smax, sew), 1
    if r < smin:
        return U(smin, sew), 1
    return U(r, sew), 0

# ===================================================================
# Scaling shift right with rounding
# ===================================================================

def ssrl(a: int, b: int, sew: int, vxrm: int = 0) -> int:
    shift = b & (sew - 1)
    return U(roundoff_unsigned(U(a, sew), shift, vxrm), sew)

def ssra(a: int, b: int, sew: int, vxrm: int = 0) -> int:
    shift = b & (sew - 1)
    return U(roundoff_signed(S(a, sew), shift, vxrm), sew)

# ===================================================================
# Narrowing clip with rounding and saturation  →  (result, vxsat)
#   a is 2*SEW, shift by b, clip to SEW
# ===================================================================

def nclipu(a: int, b: int, sew: int, vxrm: int = 0) -> tuple[int, int]:
    shift = b & (2 * sew - 1)
    r = roundoff_unsigned(U(a, 2 * sew), shift, vxrm)
    if r > M(sew):
        return M(sew), 1
    return U(r, sew), 0

def nclip(a: int, b: int, sew: int, vxrm: int = 0) -> tuple[int, int]:
    shift = b & (2 * sew - 1)
    r = roundoff_signed(S(a, 2 * sew), shift, vxrm)
    smax = (1 << (sew - 1)) - 1
    smin = -(1 << (sew - 1))
    if r > smax:
        return U(smax, sew), 1
    if r < smin:
        return U(smin, sew), 1
    return U(r, sew), 0
