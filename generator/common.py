from __future__ import annotations

"""Common constants and utility functions for the RVV test generator."""

import struct

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SEWS: tuple[int, ...] = (8, 16, 32, 64)
"""All SEW values defined in RVV 1.0."""

FP_SEWS: tuple[int, ...] = (32, 64)
"""FP-capable SEW values for standard V extension (no Zvfh)."""

WIDENING_SEWS: tuple[int, ...] = (8, 16, 32)
"""Source SEWs valid for widening operations (dest SEW ≤ 64)."""

NUM_ELEMS: int = 4
"""Number of test elements per test case.  Fits any VLEN ≥ 128 with LMUL=1."""

DEFAULT_LMUL: str = "m1"

# ---------------------------------------------------------------------------
# Bit-manipulation helpers
# ---------------------------------------------------------------------------

def M(sew: int) -> int:
    """Bitmask for an unsigned *sew*-bit value."""
    return (1 << sew) - 1


def S(val: int, sew: int) -> int:
    """Sign-extend an *sew*-bit value to a Python (arbitrary-width) int."""
    val = val & M(sew)
    if val & (1 << (sew - 1)):
        return val - (1 << sew)
    return val


def U(val: int, sew: int) -> int:
    """Truncate *val* to an unsigned *sew*-bit value."""
    return val & M(sew)

# ---------------------------------------------------------------------------
# Fixed-point rounding helper (used by several fixed-point instructions)
# ---------------------------------------------------------------------------

def roundoff_unsigned(v: int, d: int, vxrm: int) -> int:
    """Unsigned rounding: (v >> d) + r, per RVV spec Table 5."""
    if d == 0:
        return v
    if vxrm == 0:       # rnu – round-to-nearest-up
        r = (v >> (d - 1)) & 1
    elif vxrm == 1:     # rne – round-to-nearest-even
        lo = v & ((1 << (d - 1)) - 1)
        r = ((v >> (d - 1)) & 1) & ((lo != 0) | ((v >> d) & 1))
    elif vxrm == 2:     # rdn – truncate
        r = 0
    else:               # rod – round-to-odd (jam)
        lo = v & ((1 << d) - 1)
        r = int(not ((v >> d) & 1)) & int(lo != 0)
    return (v >> d) + r


def roundoff_signed(v: int, d: int, vxrm: int) -> int:
    """Signed rounding: arithmetic shift + rounding increment."""
    if d == 0:
        return v
    # Extract the needed bits from *v* directly (works for negative v too
    # because Python ints have arbitrary precision and >> is arithmetic).
    if vxrm == 0:
        r = (v >> (d - 1)) & 1
    elif vxrm == 1:
        lo = v & ((1 << (d - 1)) - 1)
        r = ((v >> (d - 1)) & 1) & ((lo != 0) | ((v >> d) & 1))
    elif vxrm == 2:
        r = 0
    else:
        lo = v & ((1 << d) - 1)
        r = int(not ((v >> d) & 1)) & int(lo != 0)
    return (v >> d) + r

# ---------------------------------------------------------------------------
# IEEE 754 helpers
# ---------------------------------------------------------------------------

def f32_to_bits(f: float) -> int:
    """Convert Python float → IEEE 754 binary32 bit-pattern (uint32)."""
    return struct.unpack("<I", struct.pack("<f", f))[0]


def bits_to_f32(i: int) -> float:
    """Convert uint32 bit-pattern → Python float (via binary32)."""
    return struct.unpack("<f", struct.pack("<I", i & 0xFFFF_FFFF))[0]


def f64_to_bits(f: float) -> int:
    """Convert Python float → IEEE 754 binary64 bit-pattern (uint64)."""
    return struct.unpack("<Q", struct.pack("<d", f))[0]


def bits_to_f64(i: int) -> float:
    """Convert uint64 bit-pattern → Python float (via binary64)."""
    return struct.unpack("<d", struct.pack("<Q", i & 0xFFFF_FFFF_FFFF_FFFF))[0]

# ---------------------------------------------------------------------------
# Assembly-level formatting helpers
# ---------------------------------------------------------------------------

_DATA_DIRECTIVES: dict[int, str] = {8: ".byte", 16: ".half", 32: ".word", 64: ".dword"}

def data_directive(sew: int) -> str:
    return _DATA_DIRECTIVES[sew]


def format_hex(val: int, sew: int) -> str:
    """Format *val* as a ``0x…`` hex literal of appropriate width."""
    return f"0x{val & M(sew):0{sew // 4}x}"


def format_data_line(values: list[int], sew: int) -> str:
    """Return e.g. ``    .word 0x00000001, 0x00000002, …``."""
    directive = data_directive(sew)
    hexvals = ", ".join(format_hex(v, sew) for v in values)
    return f"    {directive} {hexvals}"


def witness_pattern(sew: int) -> list[int]:
    """Distinctive per-SEW patterns used for witness registers."""
    return {
        8:  [0xDE, 0xAD, 0xBE, 0xEF],
        16: [0xDEAD, 0xBEEF, 0xCAFE, 0xBABE],
        32: [0xDEAD_BEEF, 0xCAFE_BABE, 0x1234_5678, 0x9ABC_DEF0],
        64: [0xDEAD_BEEF_CAFE_BABE, 0x1234_5678_9ABC_DEF0,
             0xFEDC_BA98_7654_3210, 0x0123_4567_89AB_CDEF],
    }[sew]
