/* Auto-generated test for vssubu.vv
 * Saturating vssubu.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vssubu.vv e8 basic: wrong result
 *     2 = vssubu.vv e8 basic: witness changed
 *     3 = vssubu.vv e8 basic: CSR side-effect
 *     4 = vssubu.vv e8 zeros: wrong result
 *     5 = vssubu.vv e8 zeros: witness changed
 *     6 = vssubu.vv e8 zeros: CSR side-effect
 *     7 = vssubu.vv e8 max: wrong result
 *     8 = vssubu.vv e8 max: witness changed
 *     9 = vssubu.vv e8 max: CSR side-effect
 *    10 = vssubu.vv e8 signed: wrong result
 *    11 = vssubu.vv e8 signed: witness changed
 *    12 = vssubu.vv e8 signed: CSR side-effect
 *    13 = vssubu.vv e8 mixed: wrong result
 *    14 = vssubu.vv e8 mixed: witness changed
 *    15 = vssubu.vv e8 mixed: CSR side-effect
 *    16 = vssubu.vv e8 extreme: wrong result
 *    17 = vssubu.vv e8 extreme: witness changed
 *    18 = vssubu.vv e8 extreme: CSR side-effect
 *    19 = vssubu.vv e16 basic: wrong result
 *    20 = vssubu.vv e16 basic: witness changed
 *    21 = vssubu.vv e16 basic: CSR side-effect
 *    22 = vssubu.vv e16 zeros: wrong result
 *    23 = vssubu.vv e16 zeros: witness changed
 *    24 = vssubu.vv e16 zeros: CSR side-effect
 *    25 = vssubu.vv e16 max: wrong result
 *    26 = vssubu.vv e16 max: witness changed
 *    27 = vssubu.vv e16 max: CSR side-effect
 *    28 = vssubu.vv e16 signed: wrong result
 *    29 = vssubu.vv e16 signed: witness changed
 *    30 = vssubu.vv e16 signed: CSR side-effect
 *    31 = vssubu.vv e16 mixed: wrong result
 *    32 = vssubu.vv e16 mixed: witness changed
 *    33 = vssubu.vv e16 mixed: CSR side-effect
 *    34 = vssubu.vv e16 extreme: wrong result
 *    35 = vssubu.vv e16 extreme: witness changed
 *    36 = vssubu.vv e16 extreme: CSR side-effect
 *    37 = vssubu.vv e32 basic: wrong result
 *    38 = vssubu.vv e32 basic: witness changed
 *    39 = vssubu.vv e32 basic: CSR side-effect
 *    40 = vssubu.vv e32 zeros: wrong result
 *    41 = vssubu.vv e32 zeros: witness changed
 *    42 = vssubu.vv e32 zeros: CSR side-effect
 *    43 = vssubu.vv e32 max: wrong result
 *    44 = vssubu.vv e32 max: witness changed
 *    45 = vssubu.vv e32 max: CSR side-effect
 *    46 = vssubu.vv e32 signed: wrong result
 *    47 = vssubu.vv e32 signed: witness changed
 *    48 = vssubu.vv e32 signed: CSR side-effect
 *    49 = vssubu.vv e32 mixed: wrong result
 *    50 = vssubu.vv e32 mixed: witness changed
 *    51 = vssubu.vv e32 mixed: CSR side-effect
 *    52 = vssubu.vv e32 extreme: wrong result
 *    53 = vssubu.vv e32 extreme: witness changed
 *    54 = vssubu.vv e32 extreme: CSR side-effect
 *    55 = vssubu.vv e32 basic overlap(vd=vs2): wrong result
 *    56 = vssubu.vv e32 basic overlap(vd=vs2): witness changed
 *    57 = vssubu.vv e32 basic overlap(vd=vs2): CSR side-effect
 *    58 = vssubu.vv e32 basic masked: wrong result
 *    59 = vssubu.vv e32 basic masked: witness changed
 *    60 = vssubu.vv e32 basic masked: CSR side-effect
 *    61 = vssubu.vv e64 basic: wrong result
 *    62 = vssubu.vv e64 basic: witness changed
 *    63 = vssubu.vv e64 basic: CSR side-effect
 *    64 = vssubu.vv e64 zeros: wrong result
 *    65 = vssubu.vv e64 zeros: witness changed
 *    66 = vssubu.vv e64 zeros: CSR side-effect
 *    67 = vssubu.vv e64 max: wrong result
 *    68 = vssubu.vv e64 max: witness changed
 *    69 = vssubu.vv e64 max: CSR side-effect
 *    70 = vssubu.vv e64 signed: wrong result
 *    71 = vssubu.vv e64 signed: witness changed
 *    72 = vssubu.vv e64 signed: CSR side-effect
 *    73 = vssubu.vv e64 mixed: wrong result
 *    74 = vssubu.vv e64 mixed: witness changed
 *    75 = vssubu.vv e64 mixed: CSR side-effect
 *    76 = vssubu.vv e64 extreme: wrong result
 *    77 = vssubu.vv e64 extreme: witness changed
 *    78 = vssubu.vv e64 extreme: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vssubu.vv SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    la t1, tc1_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    SET_TEST_NUM 2
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 4
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 4-6: vssubu.vv SEW=8 zeros */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    la t1, tc4_s1
    vle8.v v20, (t1)
    la t1, tc4_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    SET_TEST_NUM 5
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 4
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 7-9: vssubu.vv SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_s2
    vle8.v v16, (t1)
    la t1, tc7_s1
    vle8.v v20, (t1)
    la t1, tc7_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 4
    SET_TEST_NUM 8
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 4
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 10-12: vssubu.vv SEW=8 signed */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc10_s2
    vle8.v v16, (t1)
    la t1, tc10_s1
    vle8.v v20, (t1)
    la t1, tc10_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 10
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 4
    SET_TEST_NUM 11
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 4
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 13-15: vssubu.vv SEW=8 mixed */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc13_s2
    vle8.v v16, (t1)
    la t1, tc13_s1
    vle8.v v20, (t1)
    la t1, tc13_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 13
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 4
    SET_TEST_NUM 14
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 4
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 16-18: vssubu.vv SEW=8 extreme */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc16_s2
    vle8.v v16, (t1)
    la t1, tc16_s1
    vle8.v v20, (t1)
    la t1, tc16_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 16
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 4
    SET_TEST_NUM 17
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc16_w, 4
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 19-21: vssubu.vv SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc19_s2
    vle16.v v16, (t1)
    la t1, tc19_s1
    vle16.v v20, (t1)
    la t1, tc19_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 19
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 8
    SET_TEST_NUM 20
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc19_w, 8
    SET_TEST_NUM 21
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 22-24: vssubu.vv SEW=16 zeros */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc22_s2
    vle16.v v16, (t1)
    la t1, tc22_s1
    vle16.v v20, (t1)
    la t1, tc22_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 22
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc22_exp, 8
    SET_TEST_NUM 23
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc22_w, 8
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 25-27: vssubu.vv SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc25_s2
    vle16.v v16, (t1)
    la t1, tc25_s1
    vle16.v v20, (t1)
    la t1, tc25_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 25
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc25_exp, 8
    SET_TEST_NUM 26
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc25_w, 8
    SET_TEST_NUM 27
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 28-30: vssubu.vv SEW=16 signed */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc28_s2
    vle16.v v16, (t1)
    la t1, tc28_s1
    vle16.v v20, (t1)
    la t1, tc28_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 28
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc28_exp, 8
    SET_TEST_NUM 29
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc28_w, 8
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 31-33: vssubu.vv SEW=16 mixed */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc31_s2
    vle16.v v16, (t1)
    la t1, tc31_s1
    vle16.v v20, (t1)
    la t1, tc31_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 31
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc31_exp, 8
    SET_TEST_NUM 32
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc31_w, 8
    SET_TEST_NUM 33
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 34-36: vssubu.vv SEW=16 extreme */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc34_s2
    vle16.v v16, (t1)
    la t1, tc34_s1
    vle16.v v20, (t1)
    la t1, tc34_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 34
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc34_exp, 8
    SET_TEST_NUM 35
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc34_w, 8
    SET_TEST_NUM 36
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 37-39: vssubu.vv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc37_s2
    vle32.v v16, (t1)
    la t1, tc37_s1
    vle32.v v20, (t1)
    la t1, tc37_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 37
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc37_exp, 16
    SET_TEST_NUM 38
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc37_w, 16
    SET_TEST_NUM 39
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 40-42: vssubu.vv SEW=32 zeros */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc40_s2
    vle32.v v16, (t1)
    la t1, tc40_s1
    vle32.v v20, (t1)
    la t1, tc40_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 40
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc40_exp, 16
    SET_TEST_NUM 41
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc40_w, 16
    SET_TEST_NUM 42
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 43-45: vssubu.vv SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc43_s2
    vle32.v v16, (t1)
    la t1, tc43_s1
    vle32.v v20, (t1)
    la t1, tc43_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 43
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc43_exp, 16
    SET_TEST_NUM 44
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc43_w, 16
    SET_TEST_NUM 45
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 46-48: vssubu.vv SEW=32 signed */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc46_s2
    vle32.v v16, (t1)
    la t1, tc46_s1
    vle32.v v20, (t1)
    la t1, tc46_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 46
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc46_exp, 16
    SET_TEST_NUM 47
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc46_w, 16
    SET_TEST_NUM 48
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 49-51: vssubu.vv SEW=32 mixed */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc49_s2
    vle32.v v16, (t1)
    la t1, tc49_s1
    vle32.v v20, (t1)
    la t1, tc49_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 49
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc49_exp, 16
    SET_TEST_NUM 50
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc49_w, 16
    SET_TEST_NUM 51
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 52-54: vssubu.vv SEW=32 extreme */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc52_s2
    vle32.v v16, (t1)
    la t1, tc52_s1
    vle32.v v20, (t1)
    la t1, tc52_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 52
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc52_exp, 16
    SET_TEST_NUM 53
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc52_w, 16
    SET_TEST_NUM 54
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 55-57: vssubu.vv SEW=32 basic (vd=vs2) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc55_s2
    vle32.v v8, (t1)
    la t1, tc55_s1
    vle32.v v20, (t1)
    la t1, tc55_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v8, v20
    SET_TEST_NUM 55
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc55_exp, 16
    SET_TEST_NUM 56
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc55_w, 16
    SET_TEST_NUM 57
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 58-60: vssubu.vv SEW=32 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc58_s2
    vle32.v v16, (t1)
    la t1, tc58_s1
    vle32.v v20, (t1)
    la t1, tc58_vd
    vle32.v v8, (t1)
    la t1, tc58_mask
    vlm.v v0, (t1)
    la t1, tc58_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20, v0.t
    SET_TEST_NUM 58
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc58_exp, 16
    SET_TEST_NUM 59
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc58_w, 16
    SET_TEST_NUM 60
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 61-63: vssubu.vv SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc61_s2
    vle64.v v16, (t1)
    la t1, tc61_s1
    vle64.v v20, (t1)
    la t1, tc61_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 61
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc61_exp, 32
    SET_TEST_NUM 62
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc61_w, 32
    SET_TEST_NUM 63
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 64-66: vssubu.vv SEW=64 zeros */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc64_s2
    vle64.v v16, (t1)
    la t1, tc64_s1
    vle64.v v20, (t1)
    la t1, tc64_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 64
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc64_exp, 32
    SET_TEST_NUM 65
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc64_w, 32
    SET_TEST_NUM 66
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 67-69: vssubu.vv SEW=64 max */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc67_s2
    vle64.v v16, (t1)
    la t1, tc67_s1
    vle64.v v20, (t1)
    la t1, tc67_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 67
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc67_exp, 32
    SET_TEST_NUM 68
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc67_w, 32
    SET_TEST_NUM 69
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 70-72: vssubu.vv SEW=64 signed */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc70_s2
    vle64.v v16, (t1)
    la t1, tc70_s1
    vle64.v v20, (t1)
    la t1, tc70_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 70
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc70_exp, 32
    SET_TEST_NUM 71
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc70_w, 32
    SET_TEST_NUM 72
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 73-75: vssubu.vv SEW=64 mixed */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc73_s2
    vle64.v v16, (t1)
    la t1, tc73_s1
    vle64.v v20, (t1)
    la t1, tc73_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 73
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc73_exp, 32
    SET_TEST_NUM 74
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc73_w, 32
    SET_TEST_NUM 75
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Test 76-78: vssubu.vv SEW=64 extreme */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc76_s2
    vle64.v v16, (t1)
    la t1, tc76_s1
    vle64.v v20, (t1)
    la t1, tc76_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vssubu.vv v8, v16, v20
    SET_TEST_NUM 76
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc76_exp, 32
    SET_TEST_NUM 77
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc76_w, 32
    SET_TEST_NUM 78
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_s1:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .byte 0x00, 0x00, 0x00, 0x00
tc1_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc4_s2:
    .byte 0x00, 0x00, 0x00, 0x00
tc4_s1:
    .byte 0x00, 0x00, 0x00, 0x00
tc4_exp:
    .byte 0x00, 0x00, 0x00, 0x00
tc4_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc7_s2:
    .byte 0xff, 0xff, 0xff, 0xff
tc7_s1:
    .byte 0x01, 0x02, 0x03, 0x04
tc7_exp:
    .byte 0xfe, 0xfd, 0xfc, 0xfb
tc7_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc10_s2:
    .byte 0x7f, 0x80, 0x81, 0xff
tc10_s1:
    .byte 0x01, 0x01, 0x01, 0x01
tc10_exp:
    .byte 0x7e, 0x7f, 0x80, 0xfe
tc10_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc13_s2:
    .byte 0x00, 0x01, 0x80, 0xff
tc13_s1:
    .byte 0xff, 0x80, 0x01, 0x00
tc13_exp:
    .byte 0x00, 0x00, 0x7f, 0xff
tc13_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc16_s2:
    .byte 0x80, 0xff, 0x80, 0xff
tc16_s1:
    .byte 0x80, 0xff, 0x00, 0xff
tc16_exp:
    .byte 0x00, 0x00, 0x80, 0x00
tc16_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc19_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc19_s1:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc19_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc19_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc22_s2:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc22_s1:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc22_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc22_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc25_s2:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc25_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc25_exp:
    .half 0xfffe, 0xfffd, 0xfffc, 0xfffb
tc25_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc28_s2:
    .half 0x7fff, 0x8000, 0x8001, 0xffff
tc28_s1:
    .half 0x0001, 0x0001, 0x0001, 0x0001
tc28_exp:
    .half 0x7ffe, 0x7fff, 0x8000, 0xfffe
tc28_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc31_s2:
    .half 0x0000, 0x0001, 0x8000, 0xffff
tc31_s1:
    .half 0xffff, 0x8000, 0x0001, 0x0000
tc31_exp:
    .half 0x0000, 0x0000, 0x7fff, 0xffff
tc31_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc34_s2:
    .half 0x8000, 0xffff, 0x8000, 0xffff
tc34_s1:
    .half 0x8000, 0xffff, 0x0000, 0xffff
tc34_exp:
    .half 0x0000, 0x0000, 0x8000, 0x0000
tc34_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 2
tc37_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc37_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc37_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc37_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc40_s2:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc40_s1:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc40_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc40_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc43_s2:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc43_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc43_exp:
    .word 0xfffffffe, 0xfffffffd, 0xfffffffc, 0xfffffffb
tc43_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc46_s2:
    .word 0x7fffffff, 0x80000000, 0x80000001, 0xffffffff
tc46_s1:
    .word 0x00000001, 0x00000001, 0x00000001, 0x00000001
tc46_exp:
    .word 0x7ffffffe, 0x7fffffff, 0x80000000, 0xfffffffe
tc46_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc49_s2:
    .word 0x00000000, 0x00000001, 0x80000000, 0xffffffff
tc49_s1:
    .word 0xffffffff, 0x80000000, 0x00000001, 0x00000000
tc49_exp:
    .word 0x00000000, 0x00000000, 0x7fffffff, 0xffffffff
tc49_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc52_s2:
    .word 0x80000000, 0xffffffff, 0x80000000, 0xffffffff
tc52_s1:
    .word 0x80000000, 0xffffffff, 0x00000000, 0xffffffff
tc52_exp:
    .word 0x00000000, 0x00000000, 0x80000000, 0x00000000
tc52_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc55_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc55_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc55_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc55_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 1
tc58_mask:
    .byte 10
.align 2
tc58_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc58_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc58_vd:
    .word 0x0000dead, 0x0000dead, 0x0000dead, 0x0000dead
tc58_exp:
    .word 0x0000dead, 0x00000000, 0x0000dead, 0x00000000
tc58_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc61_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc61_s1:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc61_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc61_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc64_s2:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc64_s1:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc64_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc64_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc67_s2:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
tc67_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc67_exp:
    .dword 0xfffffffffffffffe, 0xfffffffffffffffd, 0xfffffffffffffffc, 0xfffffffffffffffb
tc67_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc70_s2:
    .dword 0x7fffffffffffffff, 0x8000000000000000, 0x8000000000000001, 0xffffffffffffffff
tc70_s1:
    .dword 0x0000000000000001, 0x0000000000000001, 0x0000000000000001, 0x0000000000000001
tc70_exp:
    .dword 0x7ffffffffffffffe, 0x7fffffffffffffff, 0x8000000000000000, 0xfffffffffffffffe
tc70_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc73_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0x8000000000000000, 0xffffffffffffffff
tc73_s1:
    .dword 0xffffffffffffffff, 0x8000000000000000, 0x0000000000000001, 0x0000000000000000
tc73_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x7fffffffffffffff, 0xffffffffffffffff
tc73_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc76_s2:
    .dword 0x8000000000000000, 0xffffffffffffffff, 0x8000000000000000, 0xffffffffffffffff
tc76_s1:
    .dword 0x8000000000000000, 0xffffffffffffffff, 0x0000000000000000, 0xffffffffffffffff
tc76_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 0x0000000000000000
tc76_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

