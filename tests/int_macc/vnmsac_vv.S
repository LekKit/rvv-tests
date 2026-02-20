/* Auto-generated test for vnmsac.vv
 * Multiply-accumulate vnmsac.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vnmsac.vv e8 basic: wrong result
 *     2 = vnmsac.vv e8 basic: CSR side-effect
 *     3 = vnmsac.vv e8 accum: wrong result
 *     4 = vnmsac.vv e8 accum: CSR side-effect
 *     5 = vnmsac.vv e8 max: wrong result
 *     6 = vnmsac.vv e8 max: CSR side-effect
 *     7 = vnmsac.vv e16 basic: wrong result
 *     8 = vnmsac.vv e16 basic: CSR side-effect
 *     9 = vnmsac.vv e16 accum: wrong result
 *    10 = vnmsac.vv e16 accum: CSR side-effect
 *    11 = vnmsac.vv e16 max: wrong result
 *    12 = vnmsac.vv e16 max: CSR side-effect
 *    13 = vnmsac.vv e32 basic: wrong result
 *    14 = vnmsac.vv e32 basic: CSR side-effect
 *    15 = vnmsac.vv e32 accum: wrong result
 *    16 = vnmsac.vv e32 accum: CSR side-effect
 *    17 = vnmsac.vv e32 max: wrong result
 *    18 = vnmsac.vv e32 max: CSR side-effect
 *    19 = vnmsac.vv e64 basic: wrong result
 *    20 = vnmsac.vv e64 basic: CSR side-effect
 *    21 = vnmsac.vv e64 accum: wrong result
 *    22 = vnmsac.vv e64 accum: CSR side-effect
 *    23 = vnmsac.vv e64 max: wrong result
 *    24 = vnmsac.vv e64 max: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vnmsac.vv SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_vd
    vle8.v v8, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    la t1, tc1_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    /* Test 3-4: vnmsac.vv SEW=8 accum */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_vd
    vle8.v v8, (t1)
    la t1, tc3_s1
    vle8.v v20, (t1)
    la t1, tc3_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED

    /* Test 5-6: vnmsac.vv SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_vd
    vle8.v v8, (t1)
    la t1, tc5_s1
    vle8.v v20, (t1)
    la t1, tc5_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 5
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 4
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-8: vnmsac.vv SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_vd
    vle16.v v8, (t1)
    la t1, tc7_s1
    vle16.v v20, (t1)
    la t1, tc7_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 7
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED

    /* Test 9-10: vnmsac.vv SEW=16 accum */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_vd
    vle16.v v8, (t1)
    la t1, tc9_s1
    vle16.v v20, (t1)
    la t1, tc9_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 9
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 8
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED

    /* Test 11-12: vnmsac.vv SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_vd
    vle16.v v8, (t1)
    la t1, tc11_s1
    vle16.v v20, (t1)
    la t1, tc11_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 11
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 8
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13-14: vnmsac.vv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_vd
    vle32.v v8, (t1)
    la t1, tc13_s1
    vle32.v v20, (t1)
    la t1, tc13_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16
    SET_TEST_NUM 14
    CHECK_CSRS_UNCHANGED

    /* Test 15-16: vnmsac.vv SEW=32 accum */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc15_vd
    vle32.v v8, (t1)
    la t1, tc15_s1
    vle32.v v20, (t1)
    la t1, tc15_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 15
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 16
    SET_TEST_NUM 16
    CHECK_CSRS_UNCHANGED

    /* Test 17-18: vnmsac.vv SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc17_vd
    vle32.v v8, (t1)
    la t1, tc17_s1
    vle32.v v20, (t1)
    la t1, tc17_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 17
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 16
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    /* Test 19-20: vnmsac.vv SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc19_vd
    vle64.v v8, (t1)
    la t1, tc19_s1
    vle64.v v20, (t1)
    la t1, tc19_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 19
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 32
    SET_TEST_NUM 20
    CHECK_CSRS_UNCHANGED

    /* Test 21-22: vnmsac.vv SEW=64 accum */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc21_vd
    vle64.v v8, (t1)
    la t1, tc21_s1
    vle64.v v20, (t1)
    la t1, tc21_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 21
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc21_exp, 32
    SET_TEST_NUM 22
    CHECK_CSRS_UNCHANGED

    /* Test 23-24: vnmsac.vv SEW=64 max */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc23_vd
    vle64.v v8, (t1)
    la t1, tc23_s1
    vle64.v v20, (t1)
    la t1, tc23_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vnmsac.vv v8, v20, v16
    SET_TEST_NUM 23
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc23_exp, 32
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_vd:
    .byte 0x00, 0x00, 0x00, 0x00
tc1_s1:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_s2:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .byte 0xf6, 0xd8, 0xa6, 0x60
.align 1
tc3_vd:
    .byte 0x64, 0xc8, 0x2c, 0x90
tc3_s1:
    .byte 0x01, 0x01, 0x01, 0x01
tc3_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc3_exp:
    .byte 0x63, 0xc6, 0x29, 0x8c
.align 1
tc5_vd:
    .byte 0x00, 0x00, 0x00, 0x00
tc5_s1:
    .byte 0xff, 0xff, 0x01, 0x00
tc5_s2:
    .byte 0x01, 0x01, 0xff, 0xff
tc5_exp:
    .byte 0x01, 0x01, 0x01, 0x00
.align 1
tc7_vd:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc7_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc7_s2:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc7_exp:
    .half 0xfff6, 0xffd8, 0xffa6, 0xff60
.align 1
tc9_vd:
    .half 0x0064, 0x00c8, 0x012c, 0x0190
tc9_s1:
    .half 0x0001, 0x0001, 0x0001, 0x0001
tc9_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc9_exp:
    .half 0x0063, 0x00c6, 0x0129, 0x018c
.align 1
tc11_vd:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc11_s1:
    .half 0xffff, 0xffff, 0x0001, 0x0000
tc11_s2:
    .half 0x0001, 0x0001, 0xffff, 0xffff
tc11_exp:
    .half 0x0001, 0x0001, 0x0001, 0x0000
.align 2
tc13_vd:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc13_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc13_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc13_exp:
    .word 0xfffffff6, 0xffffffd8, 0xffffffa6, 0xffffff60
.align 2
tc15_vd:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
tc15_s1:
    .word 0x00000001, 0x00000001, 0x00000001, 0x00000001
tc15_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc15_exp:
    .word 0x00000063, 0x000000c6, 0x00000129, 0x0000018c
.align 2
tc17_vd:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc17_s1:
    .word 0xffffffff, 0xffffffff, 0x00000001, 0x00000000
tc17_s2:
    .word 0x00000001, 0x00000001, 0xffffffff, 0xffffffff
tc17_exp:
    .word 0x00000001, 0x00000001, 0x00000001, 0x00000000
.align 3
tc19_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc19_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc19_s2:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc19_exp:
    .dword 0xfffffffffffffff6, 0xffffffffffffffd8, 0xffffffffffffffa6, 0xffffffffffffff60
.align 3
tc21_vd:
    .dword 0x0000000000000064, 0x00000000000000c8, 0x000000000000012c, 0x0000000000000190
tc21_s1:
    .dword 0x0000000000000001, 0x0000000000000001, 0x0000000000000001, 0x0000000000000001
tc21_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc21_exp:
    .dword 0x0000000000000063, 0x00000000000000c6, 0x0000000000000129, 0x000000000000018c
.align 3
tc23_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc23_s1:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0x0000000000000001, 0x0000000000000000
tc23_s2:
    .dword 0x0000000000000001, 0x0000000000000001, 0xffffffffffffffff, 0xffffffffffffffff
tc23_exp:
    .dword 0x0000000000000001, 0x0000000000000001, 0x0000000000000001, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

