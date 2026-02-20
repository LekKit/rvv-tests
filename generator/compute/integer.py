from __future__ import annotations

"""Compute functions for all integer vector instructions."""

import math

from ..common import M, S, U

# ===================================================================
# Single-width arithmetic
# ===================================================================

def add(a: int, b: int, sew: int) -> int:
    return U(a + b, sew)

def sub(a: int, b: int, sew: int) -> int:
    return U(a - b, sew)

def rsub(a: int, b: int, sew: int) -> int:
    return U(b - a, sew)

# ===================================================================
# Bitwise logical
# ===================================================================

def and_(a: int, b: int, sew: int) -> int:
    return (a & b) & M(sew)

def or_(a: int, b: int, sew: int) -> int:
    return (a | b) & M(sew)

def xor(a: int, b: int, sew: int) -> int:
    return (a ^ b) & M(sew)

# ===================================================================
# Shifts
# ===================================================================

def sll(a: int, b: int, sew: int) -> int:
    return U(a << (b & (sew - 1)), sew)

def srl(a: int, b: int, sew: int) -> int:
    return U(a, sew) >> (b & (sew - 1))

def sra(a: int, b: int, sew: int) -> int:
    return U(S(a, sew) >> (b & (sew - 1)), sew)

# ===================================================================
# Min / Max
# ===================================================================

def minu(a: int, b: int, sew: int) -> int:
    return min(U(a, sew), U(b, sew))

def min_(a: int, b: int, sew: int) -> int:
    return U(min(S(a, sew), S(b, sew)), sew)

def maxu(a: int, b: int, sew: int) -> int:
    return max(U(a, sew), U(b, sew))

def max_(a: int, b: int, sew: int) -> int:
    return U(max(S(a, sew), S(b, sew)), sew)

# ===================================================================
# Multiply
# ===================================================================

def mul(a: int, b: int, sew: int) -> int:
    return U(S(a, sew) * S(b, sew), sew)

def mulh(a: int, b: int, sew: int) -> int:
    return U((S(a, sew) * S(b, sew)) >> sew, sew)

def mulhu(a: int, b: int, sew: int) -> int:
    return U((U(a, sew) * U(b, sew)) >> sew, sew)

def mulhsu(a: int, b: int, sew: int) -> int:
    """vs2 is signed, vs1/rs1 is unsigned."""
    return U((S(a, sew) * U(b, sew)) >> sew, sew)

# ===================================================================
# Divide / Remainder
# ===================================================================

def _trunc_div(a: int, b: int) -> int:
    """Truncating integer division (toward zero), matching C / RISC-V."""
    if b == 0:
        raise ZeroDivisionError
    # Python's // floors; we need truncation.
    return int(math.copysign(1, a * b)) * (abs(a) // abs(b)) if a * b else 0

def divu(a: int, b: int, sew: int) -> int:
    ua, ub = U(a, sew), U(b, sew)
    return M(sew) if ub == 0 else U(ua // ub, sew)

def div(a: int, b: int, sew: int) -> int:
    sa, sb = S(a, sew), S(b, sew)
    if sb == 0:
        return M(sew)
    smin = -(1 << (sew - 1))
    if sa == smin and sb == -1:
        return U(smin, sew)  # overflow → same bit pattern
    q = -(abs(sa) // abs(sb)) if (sa < 0) != (sb < 0) else abs(sa) // abs(sb)
    return U(q, sew)

def remu(a: int, b: int, sew: int) -> int:
    ua, ub = U(a, sew), U(b, sew)
    return ua if ub == 0 else U(ua % ub, sew)

def rem(a: int, b: int, sew: int) -> int:
    sa, sb = S(a, sew), S(b, sew)
    if sb == 0:
        return U(sa, sew)
    smin = -(1 << (sew - 1))
    if sa == smin and sb == -1:
        return 0
    q = -(abs(sa) // abs(sb)) if (sa < 0) != (sb < 0) else abs(sa) // abs(sb)
    return U(sa - q * sb, sew)

# ===================================================================
# Compare  (return 0 or 1)
# ===================================================================

def mseq(a: int, b: int, sew: int) -> int:
    return int(U(a, sew) == U(b, sew))

def msne(a: int, b: int, sew: int) -> int:
    return int(U(a, sew) != U(b, sew))

def msltu(a: int, b: int, sew: int) -> int:
    return int(U(a, sew) < U(b, sew))

def mslt(a: int, b: int, sew: int) -> int:
    return int(S(a, sew) < S(b, sew))

def msleu(a: int, b: int, sew: int) -> int:
    return int(U(a, sew) <= U(b, sew))

def msle(a: int, b: int, sew: int) -> int:
    return int(S(a, sew) <= S(b, sew))

def msgtu(a: int, b: int, sew: int) -> int:
    return int(U(a, sew) > U(b, sew))

def msgt(a: int, b: int, sew: int) -> int:
    return int(S(a, sew) > S(b, sew))

# ===================================================================
# Widening arithmetic  (result is 2*SEW)
# ===================================================================

def waddu(a: int, b: int, sew: int) -> int:
    return U(U(a, sew) + U(b, sew), 2 * sew)

def wadd(a: int, b: int, sew: int) -> int:
    return U(S(a, sew) + S(b, sew), 2 * sew)

def wsubu(a: int, b: int, sew: int) -> int:
    return U(U(a, sew) - U(b, sew), 2 * sew)

def wsub(a: int, b: int, sew: int) -> int:
    return U(S(a, sew) - S(b, sew), 2 * sew)

# .W forms: a is 2*SEW, b is SEW
def waddu_w(a: int, b: int, sew: int) -> int:
    return U(U(a, 2 * sew) + U(b, sew), 2 * sew)

def wadd_w(a: int, b: int, sew: int) -> int:
    return U(S(a, 2 * sew) + S(b, sew), 2 * sew)

def wsubu_w(a: int, b: int, sew: int) -> int:
    return U(U(a, 2 * sew) - U(b, sew), 2 * sew)

def wsub_w(a: int, b: int, sew: int) -> int:
    return U(S(a, 2 * sew) - S(b, sew), 2 * sew)

# Widening multiply
def wmulu(a: int, b: int, sew: int) -> int:
    return U(U(a, sew) * U(b, sew), 2 * sew)

def wmul(a: int, b: int, sew: int) -> int:
    return U(S(a, sew) * S(b, sew), 2 * sew)

def wmulsu(a: int, b: int, sew: int) -> int:
    return U(S(a, sew) * U(b, sew), 2 * sew)

# ===================================================================
# Narrowing shifts  (a is 2*SEW, result is SEW)
# ===================================================================

def nsrl(a: int, b: int, sew: int) -> int:
    return U(U(a, 2 * sew) >> (b & (2 * sew - 1)), sew)

def nsra(a: int, b: int, sew: int) -> int:
    return U(S(a, 2 * sew) >> (b & (2 * sew - 1)), sew)

# ===================================================================
# Integer extension  (unary, source is narrower)
# ===================================================================

def zext2(a: int, sew: int) -> int:
    return U(a, sew // 2)

def sext2(a: int, sew: int) -> int:
    return U(S(a, sew // 2), sew)

def zext4(a: int, sew: int) -> int:
    return U(a, sew // 4)

def sext4(a: int, sew: int) -> int:
    return U(S(a, sew // 4), sew)

def zext8(a: int, sew: int) -> int:
    return U(a, sew // 8)

def sext8(a: int, sew: int) -> int:
    return U(S(a, sew // 8), sew)

# ===================================================================
# Multiply-accumulate  (vd, vs1, vs2, sew) -> result for vd
# ===================================================================

def macc(vd: int, vs1: int, vs2: int, sew: int) -> int:
    return U(S(vd, sew) + S(vs1, sew) * S(vs2, sew), sew)

def nmsac(vd: int, vs1: int, vs2: int, sew: int) -> int:
    return U(S(vd, sew) - S(vs1, sew) * S(vs2, sew), sew)

def madd(vd: int, vs1: int, vs2: int, sew: int) -> int:
    """vd = (vs1 * vd) + vs2."""
    return U(S(vs1, sew) * S(vd, sew) + S(vs2, sew), sew)

def nmsub(vd: int, vs1: int, vs2: int, sew: int) -> int:
    """vd = -(vs1 * vd) + vs2."""
    return U(-(S(vs1, sew) * S(vd, sew)) + S(vs2, sew), sew)

# Widening multiply-accumulate
def wmaccu(vd: int, vs1: int, vs2: int, sew: int) -> int:
    return U(U(vd, 2 * sew) + U(vs1, sew) * U(vs2, sew), 2 * sew)

def wmacc(vd: int, vs1: int, vs2: int, sew: int) -> int:
    return U(S(vd, 2 * sew) + S(vs1, sew) * S(vs2, sew), 2 * sew)

def wmaccsu(vd: int, vs1: int, vs2: int, sew: int) -> int:
    """signed(vs1) × unsigned(vs2) + vd."""
    return U(S(vd, 2 * sew) + S(vs1, sew) * U(vs2, sew), 2 * sew)

def wmaccus(vd: int, vs1: int, vs2: int, sew: int) -> int:
    """unsigned(rs1) × signed(vs2) + vd."""
    return U(S(vd, 2 * sew) + U(vs1, sew) * S(vs2, sew), 2 * sew)

# ===================================================================
# Mask-register logical  (operate on bit vectors)
# ===================================================================

def vmand(a: int, b: int) -> int:
    return a & b

def vmnand(a: int, b: int) -> int:
    return ~(a & b)

def vmandn(a: int, b: int) -> int:
    return a & (~b)

def vmxor(a: int, b: int) -> int:
    return a ^ b

def vmor(a: int, b: int) -> int:
    return a | b

def vmnor(a: int, b: int) -> int:
    return ~(a | b)

def vmorn(a: int, b: int) -> int:
    return a | (~b)

def vmxnor(a: int, b: int) -> int:
    return ~(a ^ b)

# ===================================================================
# Add with carry / Subtract with borrow
# ===================================================================

def adc(a: int, b: int, cin: int, sew: int) -> int:
    """vadc: vd[i] = vs2[i] + vs1[i] + carry_in."""
    return U(U(a, sew) + U(b, sew) + cin, sew)

def sbc(a: int, b: int, bin_: int, sew: int) -> int:
    """vsbc: vd[i] = vs2[i] - vs1[i] - borrow_in."""
    return U(U(a, sew) - U(b, sew) - bin_, sew)

def madc(a: int, b: int, cin: int, sew: int) -> int:
    """vmadc: carry_out of vs2 + vs1 + cin."""
    return int((U(a, sew) + U(b, sew) + cin) > M(sew))

def madc_no_carry(a: int, b: int, sew: int) -> int:
    """vmadc (no carry-in): carry_out of vs2 + vs1."""
    return int((U(a, sew) + U(b, sew)) > M(sew))

def msbc(a: int, b: int, bin_: int, sew: int) -> int:
    """vmsbc: borrow_out of vs2 - vs1 - bin."""
    return int((U(a, sew) - U(b, sew) - bin_) < 0)

def msbc_no_borrow(a: int, b: int, sew: int) -> int:
    """vmsbc (no borrow-in): borrow_out of vs2 - vs1."""
    return int(U(a, sew) < U(b, sew))
