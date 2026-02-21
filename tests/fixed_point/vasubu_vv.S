/* Auto-generated test for vasubu.vv
 * Averaging vasubu.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vasubu.vv e8 basic: result
 *     2 = vasubu.vv e8 basic: CSR
 *     3 = vasubu.vv e8 zeros: result
 *     4 = vasubu.vv e8 zeros: CSR
 *     5 = vasubu.vv e8 max: result
 *     6 = vasubu.vv e8 max: CSR
 *     7 = vasubu.vv e8 signed: result
 *     8 = vasubu.vv e8 signed: CSR
 *     9 = vasubu.vv e8 mixed: result
 *    10 = vasubu.vv e8 mixed: CSR
 *    11 = vasubu.vv e8 extreme: result
 *    12 = vasubu.vv e8 extreme: CSR
 *    13 = vasubu.vv e16 basic: result
 *    14 = vasubu.vv e16 basic: CSR
 *    15 = vasubu.vv e16 zeros: result
 *    16 = vasubu.vv e16 zeros: CSR
 *    17 = vasubu.vv e16 max: result
 *    18 = vasubu.vv e16 max: CSR
 *    19 = vasubu.vv e16 signed: result
 *    20 = vasubu.vv e16 signed: CSR
 *    21 = vasubu.vv e16 mixed: result
 *    22 = vasubu.vv e16 mixed: CSR
 *    23 = vasubu.vv e16 extreme: result
 *    24 = vasubu.vv e16 extreme: CSR
 *    25 = vasubu.vv e32 basic: result
 *    26 = vasubu.vv e32 basic: CSR
 *    27 = vasubu.vv e32 zeros: result
 *    28 = vasubu.vv e32 zeros: CSR
 *    29 = vasubu.vv e32 max: result
 *    30 = vasubu.vv e32 max: CSR
 *    31 = vasubu.vv e32 signed: result
 *    32 = vasubu.vv e32 signed: CSR
 *    33 = vasubu.vv e32 mixed: result
 *    34 = vasubu.vv e32 mixed: CSR
 *    35 = vasubu.vv e32 extreme: result
 *    36 = vasubu.vv e32 extreme: CSR
 *    37 = vasubu.vv e64 basic: result
 *    38 = vasubu.vv e64 basic: CSR
 *    39 = vasubu.vv e64 zeros: result
 *    40 = vasubu.vv e64 zeros: CSR
 *    41 = vasubu.vv e64 max: result
 *    42 = vasubu.vv e64 max: CSR
 *    43 = vasubu.vv e64 signed: result
 *    44 = vasubu.vv e64 signed: CSR
 *    45 = vasubu.vv e64 mixed: result
 *    46 = vasubu.vv e64 mixed: CSR
 *    47 = vasubu.vv e64 extreme: result
 *    48 = vasubu.vv e64 extreme: CSR
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vasubu.vv SEW=8 basic vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    /* Test 3: vasubu.vv SEW=8 zeros vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc3_s2
    vle8.v v16, (t1)
    la t1, tc3_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED

    /* Test 5: vasubu.vv SEW=8 max vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc5_s2
    vle8.v v16, (t1)
    la t1, tc5_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 5
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 4
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7: vasubu.vv SEW=8 signed vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc7_s2
    vle8.v v16, (t1)
    la t1, tc7_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 4
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED

    /* Test 9: vasubu.vv SEW=8 mixed vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc9_s2
    vle8.v v16, (t1)
    la t1, tc9_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 9
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 4
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED

    /* Test 11: vasubu.vv SEW=8 extreme vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc11_s2
    vle8.v v16, (t1)
    la t1, tc11_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 11
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 4
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13: vasubu.vv SEW=16 basic vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc13_s2
    vle16.v v16, (t1)
    la t1, tc13_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 13
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 8
    SET_TEST_NUM 14
    CHECK_CSRS_UNCHANGED

    /* Test 15: vasubu.vv SEW=16 zeros vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc15_s2
    vle16.v v16, (t1)
    la t1, tc15_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 15
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 8
    SET_TEST_NUM 16
    CHECK_CSRS_UNCHANGED

    /* Test 17: vasubu.vv SEW=16 max vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc17_s2
    vle16.v v16, (t1)
    la t1, tc17_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 17
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 8
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    /* Test 19: vasubu.vv SEW=16 signed vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc19_s2
    vle16.v v16, (t1)
    la t1, tc19_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 19
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 8
    SET_TEST_NUM 20
    CHECK_CSRS_UNCHANGED

    /* Test 21: vasubu.vv SEW=16 mixed vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc21_s2
    vle16.v v16, (t1)
    la t1, tc21_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 21
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc21_exp, 8
    SET_TEST_NUM 22
    CHECK_CSRS_UNCHANGED

    /* Test 23: vasubu.vv SEW=16 extreme vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc23_s2
    vle16.v v16, (t1)
    la t1, tc23_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 23
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc23_exp, 8
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    /* Test 25: vasubu.vv SEW=32 basic vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc25_s2
    vle32.v v16, (t1)
    la t1, tc25_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 25
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc25_exp, 16
    SET_TEST_NUM 26
    CHECK_CSRS_UNCHANGED

    /* Test 27: vasubu.vv SEW=32 zeros vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc27_s2
    vle32.v v16, (t1)
    la t1, tc27_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 27
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc27_exp, 16
    SET_TEST_NUM 28
    CHECK_CSRS_UNCHANGED

    /* Test 29: vasubu.vv SEW=32 max vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc29_s2
    vle32.v v16, (t1)
    la t1, tc29_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 29
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc29_exp, 16
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED

    /* Test 31: vasubu.vv SEW=32 signed vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc31_s2
    vle32.v v16, (t1)
    la t1, tc31_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 31
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc31_exp, 16
    SET_TEST_NUM 32
    CHECK_CSRS_UNCHANGED

    /* Test 33: vasubu.vv SEW=32 mixed vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc33_s2
    vle32.v v16, (t1)
    la t1, tc33_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 33
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc33_exp, 16
    SET_TEST_NUM 34
    CHECK_CSRS_UNCHANGED

    /* Test 35: vasubu.vv SEW=32 extreme vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc35_s2
    vle32.v v16, (t1)
    la t1, tc35_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 35
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc35_exp, 16
    SET_TEST_NUM 36
    CHECK_CSRS_UNCHANGED

    /* Test 37: vasubu.vv SEW=64 basic vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc37_s2
    vle64.v v16, (t1)
    la t1, tc37_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 37
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc37_exp, 32
    SET_TEST_NUM 38
    CHECK_CSRS_UNCHANGED

    /* Test 39: vasubu.vv SEW=64 zeros vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc39_s2
    vle64.v v16, (t1)
    la t1, tc39_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 39
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc39_exp, 32
    SET_TEST_NUM 40
    CHECK_CSRS_UNCHANGED

    /* Test 41: vasubu.vv SEW=64 max vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc41_s2
    vle64.v v16, (t1)
    la t1, tc41_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 41
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc41_exp, 32
    SET_TEST_NUM 42
    CHECK_CSRS_UNCHANGED

    /* Test 43: vasubu.vv SEW=64 signed vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc43_s2
    vle64.v v16, (t1)
    la t1, tc43_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 43
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc43_exp, 32
    SET_TEST_NUM 44
    CHECK_CSRS_UNCHANGED

    /* Test 45: vasubu.vv SEW=64 mixed vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc45_s2
    vle64.v v16, (t1)
    la t1, tc45_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 45
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc45_exp, 32
    SET_TEST_NUM 46
    CHECK_CSRS_UNCHANGED

    /* Test 47: vasubu.vv SEW=64 extreme vxrm=0 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc47_s2
    vle64.v v16, (t1)
    la t1, tc47_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vasubu.vv v8, v16, v20
    SET_TEST_NUM 47
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc47_exp, 32
    SET_TEST_NUM 48
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_s1:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .byte 0xfc, 0xf7, 0xf3, 0xee
.align 1
tc3_s2:
    .byte 0x00, 0x00, 0x00, 0x00
tc3_s1:
    .byte 0x00, 0x00, 0x00, 0x00
tc3_exp:
    .byte 0x00, 0x00, 0x00, 0x00
.align 1
tc5_s2:
    .byte 0xff, 0xff, 0xff, 0xff
tc5_s1:
    .byte 0x01, 0x02, 0x03, 0x04
tc5_exp:
    .byte 0x7f, 0x7f, 0x7e, 0x7e
.align 1
tc7_s2:
    .byte 0x7f, 0x80, 0x81, 0xff
tc7_s1:
    .byte 0x01, 0x01, 0x01, 0x01
tc7_exp:
    .byte 0x3f, 0x40, 0x40, 0x7f
.align 1
tc9_s2:
    .byte 0x00, 0x01, 0x80, 0xff
tc9_s1:
    .byte 0xff, 0x80, 0x01, 0x00
tc9_exp:
    .byte 0x81, 0xc1, 0x40, 0x80
.align 1
tc11_s2:
    .byte 0x80, 0xff, 0x80, 0xff
tc11_s1:
    .byte 0x80, 0xff, 0x00, 0xff
tc11_exp:
    .byte 0x00, 0x00, 0x40, 0x00
.align 1
tc13_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc13_s1:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc13_exp:
    .half 0xfffc, 0xfff7, 0xfff3, 0xffee
.align 1
tc15_s2:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc15_s1:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc15_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 1
tc17_s2:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc17_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc17_exp:
    .half 0x7fff, 0x7fff, 0x7ffe, 0x7ffe
.align 1
tc19_s2:
    .half 0x7fff, 0x8000, 0x8001, 0xffff
tc19_s1:
    .half 0x0001, 0x0001, 0x0001, 0x0001
tc19_exp:
    .half 0x3fff, 0x4000, 0x4000, 0x7fff
.align 1
tc21_s2:
    .half 0x0000, 0x0001, 0x8000, 0xffff
tc21_s1:
    .half 0xffff, 0x8000, 0x0001, 0x0000
tc21_exp:
    .half 0x8001, 0xc001, 0x4000, 0x8000
.align 1
tc23_s2:
    .half 0x8000, 0xffff, 0x8000, 0xffff
tc23_s1:
    .half 0x8000, 0xffff, 0x0000, 0xffff
tc23_exp:
    .half 0x0000, 0x0000, 0x4000, 0x0000
.align 2
tc25_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc25_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc25_exp:
    .word 0xfffffffc, 0xfffffff7, 0xfffffff3, 0xffffffee
.align 2
tc27_s2:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc27_s1:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc27_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
.align 2
tc29_s2:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc29_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc29_exp:
    .word 0x7fffffff, 0x7fffffff, 0x7ffffffe, 0x7ffffffe
.align 2
tc31_s2:
    .word 0x7fffffff, 0x80000000, 0x80000001, 0xffffffff
tc31_s1:
    .word 0x00000001, 0x00000001, 0x00000001, 0x00000001
tc31_exp:
    .word 0x3fffffff, 0x40000000, 0x40000000, 0x7fffffff
.align 2
tc33_s2:
    .word 0x00000000, 0x00000001, 0x80000000, 0xffffffff
tc33_s1:
    .word 0xffffffff, 0x80000000, 0x00000001, 0x00000000
tc33_exp:
    .word 0x80000001, 0xc0000001, 0x40000000, 0x80000000
.align 2
tc35_s2:
    .word 0x80000000, 0xffffffff, 0x80000000, 0xffffffff
tc35_s1:
    .word 0x80000000, 0xffffffff, 0x00000000, 0xffffffff
tc35_exp:
    .word 0x00000000, 0x00000000, 0x40000000, 0x00000000
.align 3
tc37_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc37_s1:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc37_exp:
    .dword 0xfffffffffffffffc, 0xfffffffffffffff7, 0xfffffffffffffff3, 0xffffffffffffffee
.align 3
tc39_s2:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc39_s1:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc39_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc41_s2:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
tc41_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc41_exp:
    .dword 0x7fffffffffffffff, 0x7fffffffffffffff, 0x7ffffffffffffffe, 0x7ffffffffffffffe
.align 3
tc43_s2:
    .dword 0x7fffffffffffffff, 0x8000000000000000, 0x8000000000000001, 0xffffffffffffffff
tc43_s1:
    .dword 0x0000000000000001, 0x0000000000000001, 0x0000000000000001, 0x0000000000000001
tc43_exp:
    .dword 0x3fffffffffffffff, 0x4000000000000000, 0x4000000000000000, 0x7fffffffffffffff
.align 3
tc45_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0x8000000000000000, 0xffffffffffffffff
tc45_s1:
    .dword 0xffffffffffffffff, 0x8000000000000000, 0x0000000000000001, 0x0000000000000000
tc45_exp:
    .dword 0x8000000000000001, 0xc000000000000001, 0x4000000000000000, 0x8000000000000000
.align 3
tc47_s2:
    .dword 0x8000000000000000, 0xffffffffffffffff, 0x8000000000000000, 0xffffffffffffffff
tc47_s1:
    .dword 0x8000000000000000, 0xffffffffffffffff, 0x0000000000000000, 0xffffffffffffffff
tc47_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x4000000000000000, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

