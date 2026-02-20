from __future__ import annotations

"""Test vector definitions for all instruction families.

Each function returns a list of ``(name, *operands)`` tuples.
The exact shape depends on the instruction form (VV, VX, VI, etc.).
"""

from .common import M, S, U, f32_to_bits, f64_to_bits, NUM_ELEMS

# ===================================================================
# Integer binary  (VV form):  (name, vs2[4], vs1[4])
# ===================================================================

def binop_vv(sew: int) -> list[tuple[str, list[int], list[int]]]:
    m = M(sew)
    h = 1 << (sew - 1)
    return [
        ("basic",  [1, 2, 3, 4],         [10, 20, 30, 40]),
        ("zeros",  [0, 0, 0, 0],         [0, 0, 0, 0]),
        ("max",    [m, m, m, m],         [1, 2, 3, 4]),
        ("signed", [h - 1, h, h + 1, m], [1, 1, 1, 1]),
        ("mixed",  [0, 1, h, m],         [m, h, 1, 0]),
    ]


# ===================================================================
# Integer binary  (VX form):  (name, vs2[4], scalar)
# ===================================================================

def binop_vx(sew: int) -> list[tuple[str, list[int], int]]:
    m = M(sew)
    h = 1 << (sew - 1)
    return [
        ("basic",      [1, 2, 3, 4], 10),
        ("zero_sc",    [1, 2, 3, 4], 0),
        ("max_sc",     [1, 2, 3, 4], m),
        ("max_vec",    [m, m, m, m], 1),
        ("signed_sc",  [0, 1, h, m], h),
    ]


# ===================================================================
# Integer binary  (VI form):  (name, vs2[4], imm5)
# ===================================================================

def binop_vi(sew: int, *, unsigned_imm: bool = False) -> list[tuple[str, list[int], int]]:
    m = M(sew)
    if unsigned_imm:
        # Shift-like: 0..31
        return [
            ("imm0",  [m, 1, 0, 0x55 & m], 0),
            ("imm1",  [m, 1, 0, 0x55 & m], 1),
            ("imm7",  [m, 1, 0, 0x55 & m], min(7, sew - 1)),
            ("imm15", [m, 1, 0, 0x55 & m], min(15, sew - 1)),
        ]
    # Signed: -16..15
    return [
        ("imm0",    [1, 2, 3, 4],         0),
        ("imm1",    [1, 2, 3, 4],         1),
        ("imm_n1",  [1, 2, 3, 4],        -1),
        ("imm15",   [1, 2, 3, 4],        15),
        ("imm_n16", [0, m, 1, m - 1],   -16),
    ]


# ===================================================================
# Widening  (VV form):  (name, vs2[4], vs1[4])   — source SEW
# ===================================================================

def widening_vv(sew: int) -> list[tuple[str, list[int], list[int]]]:
    m = M(sew)
    h = 1 << (sew - 1)
    return [
        ("basic",  [1, 2, 3, 4], [5, 6, 7, 8]),
        ("max",    [m, m, m, m], [m, m, m, m]),
        ("signed", [h, h + 1, 0, m], [1, m, h, 0]),
    ]


# ===================================================================
# Widening .W form:  (name, vs2_wide[4], vs1[4])
# vs2 is 2*SEW, vs1 is SEW
# ===================================================================

def widening_wv(sew: int) -> list[tuple[str, list[int], list[int]]]:
    dsew = 2 * sew
    dm = M(dsew)
    m = M(sew)
    return [
        ("basic", [1, 2, 3, 4],             [5, 6, 7, 8]),
        ("max",   [dm, dm, dm, dm],         [m, m, m, m]),
        ("mixed", [0, 1, dm >> 1, dm],      [0, 1, m, m]),
    ]


# ===================================================================
# Narrowing  (WV form):  (name, vs2_wide[4], vs1[4])
# vs2 is 2*SEW, vs1 is shift amount (SEW-width)
# ===================================================================

def narrowing_wv(sew: int) -> list[tuple[str, list[int], list[int]]]:
    dsew = 2 * sew
    dm = M(dsew)
    return [
        ("basic", [0x100, 0x200, 0x300, 0x400], [0, 0, 0, 0]),
        ("shift", [dm, dm >> 1, 1, 0],          [sew, sew // 2, 0, 1]),
        ("mixed", [dm, 0, 1 << sew, (1 << sew) - 1], [1, 0, sew - 1, 0]),
    ]


# ===================================================================
# Compare  (VV form):  same as binop but result is mask bits
# ===================================================================

def compare_vv(sew: int) -> list[tuple[str, list[int], list[int]]]:
    m = M(sew)
    h = 1 << (sew - 1)
    return [
        ("equal",    [1, 2, 3, 4],   [1, 2, 3, 4]),
        ("less",     [0, 1, 2, 3],   [1, 2, 3, 4]),
        ("greater",  [2, 3, 4, 5],   [1, 2, 3, 4]),
        ("mixed",    [0, m, h, h-1], [m, 0, h-1, h]),
    ]


# ===================================================================
# Multiply-accumulate:  (name, vd_init[4], vs1[4], vs2[4])
# ===================================================================

def macc_vv(sew: int) -> list[tuple[str, list[int], list[int], list[int]]]:
    m = M(sew)
    return [
        ("basic", [0, 0, 0, 0],     [1, 2, 3, 4],   [10, 20, 30, 40]),
        ("accum", [100, 200, 300, 400], [1, 1, 1, 1], [1, 2, 3, 4]),
        ("max",   [0, 0, 0, 0],     [m, m, 1, 0],    [1, 1, m, m]),
    ]


# ===================================================================
# Floating-point binary  (VV form):  (name, vs2[4], vs1[4])
# Values are bit patterns
# ===================================================================

def fp_binop_vv(sew: int) -> list[tuple[str, list[int], list[int]]]:
    if sew == 32:
        b = f32_to_bits
        return [
            ("basic",   [b(1.0), b(2.0), b(3.0), b(4.0)],
                        [b(0.5), b(0.25), b(0.125), b(0.0625)]),
            ("negpos",  [b(-1.0), b(-2.0), b(3.0), b(4.0)],
                        [b(1.0), b(2.0), b(-3.0), b(-4.0)]),
            ("zeros",   [b(0.0), b(-0.0), b(0.0), b(-0.0)],
                        [b(0.0), b(0.0), b(-0.0), b(-0.0)]),
        ]
    b = f64_to_bits
    return [
        ("basic",  [b(1.0), b(2.0), b(3.0), b(4.0)],
                   [b(0.5), b(0.25), b(0.125), b(0.0625)]),
        ("negpos", [b(-1.0), b(-2.0), b(3.0), b(4.0)],
                   [b(1.0), b(2.0), b(-3.0), b(-4.0)]),
    ]


# ===================================================================
# FP binary  (VF form):  (name, vs2[4], scalar_bits)
# ===================================================================

def fp_binop_vf(sew: int) -> list[tuple[str, list[int], int]]:
    if sew == 32:
        b = f32_to_bits
        return [
            ("basic",  [b(1.0), b(2.0), b(3.0), b(4.0)], b(10.0)),
            ("zero",   [b(1.0), b(2.0), b(3.0), b(4.0)], b(0.0)),
            ("neg",    [b(1.0), b(2.0), b(3.0), b(4.0)], b(-1.0)),
        ]
    b = f64_to_bits
    return [
        ("basic", [b(1.0), b(2.0), b(3.0), b(4.0)], b(10.0)),
        ("zero",  [b(1.0), b(2.0), b(3.0), b(4.0)], b(0.0)),
    ]


# ===================================================================
# FP fused multiply-add:  (name, vd_init[4], vs1[4], vs2[4])
# ===================================================================

def fp_fma_vv(sew: int) -> list[tuple[str, list[int], list[int], list[int]]]:
    if sew == 32:
        b = f32_to_bits
        return [
            ("basic", [b(0.0)] * 4,
                      [b(1.0), b(2.0), b(3.0), b(4.0)],
                      [b(10.0), b(20.0), b(30.0), b(40.0)]),
            ("accum", [b(100.0), b(200.0), b(300.0), b(400.0)],
                      [b(1.0), b(1.0), b(1.0), b(1.0)],
                      [b(1.0), b(2.0), b(3.0), b(4.0)]),
        ]
    b = f64_to_bits
    return [
        ("basic", [b(0.0)] * 4,
                  [b(1.0), b(2.0), b(3.0), b(4.0)],
                  [b(10.0), b(20.0), b(30.0), b(40.0)]),
    ]


# ===================================================================
# FP fused multiply-add:  (name, vd_init[4], scalar_bits, vs2[4])
# VF form — scalar in vs1 position
# ===================================================================

def fp_fma_vf(sew: int) -> list[tuple[str, list[int], int, list[int]]]:
    if sew == 32:
        b = f32_to_bits
        return [
            ("basic", [b(0.0)] * 4,
                      b(2.0),
                      [b(10.0), b(20.0), b(30.0), b(40.0)]),
            ("accum", [b(100.0), b(200.0), b(300.0), b(400.0)],
                      b(1.0),
                      [b(1.0), b(2.0), b(3.0), b(4.0)]),
        ]
    b = f64_to_bits
    return [
        ("basic", [b(0.0)] * 4,
                  b(2.0),
                  [b(10.0), b(20.0), b(30.0), b(40.0)]),
    ]


# ===================================================================
# Integer reductions:  (name, scalar_init, vec[4])
# ===================================================================

def reduction_int(sew: int) -> list[tuple[str, int, list[int]]]:
    m = M(sew)
    return [
        ("basic", 0, [1, 2, 3, 4]),
        ("init",  100, [1, 2, 3, 4]),
        ("max",   0, [m, m - 1, m - 2, 0]),
    ]


# ===================================================================
# Mask logical:  (name, vs2_mask_byte, vs1_mask_byte)
# For 4 elements, mask is 4 bits packed in byte
# ===================================================================

def mask_logical() -> list[tuple[str, int, int]]:
    return [
        ("all_ones",  0b1111, 0b1111),
        ("all_zeros", 0b0000, 0b0000),
        ("mixed_1",   0b1010, 0b1100),
        ("mixed_2",   0b0101, 0b0011),
        ("one_zero",  0b1111, 0b0000),
    ]


# ===================================================================
# Permutation / slide:  (name, vs2[4], offset_or_scalar)
# ===================================================================

def slide_vx(sew: int) -> list[tuple[str, list[int], int]]:
    return [
        ("by_0", [10, 20, 30, 40], 0),
        ("by_1", [10, 20, 30, 40], 1),
        ("by_2", [10, 20, 30, 40], 2),
        ("by_3", [10, 20, 30, 40], 3),
    ]

def gather_vv(sew: int) -> list[tuple[str, list[int], list[int]]]:
    m = M(sew)
    return [
        ("identity", [10, 20, 30, 40], [0, 1, 2, 3]),
        ("reverse",  [10, 20, 30, 40], [3, 2, 1, 0]),
        ("splat",    [10, 20, 30, 40], [0, 0, 0, 0]),
        ("oob",      [10, 20, 30, 40], [0, 1, 2, m]),  # out-of-bounds → 0
    ]
