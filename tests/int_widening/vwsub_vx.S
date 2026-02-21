/* Auto-generated test for vwsub.vx
 * Widening vwsub.vx
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwsub.vx e8 sc_0x01: result
 *     2 = vwsub.vx e8 sc_0x01: CSR side-effect
 *     3 = vwsub.vx e8 sc_0xff: result
 *     4 = vwsub.vx e8 sc_0xff: CSR side-effect
 *     5 = vwsub.vx e8 sc_0x00: result
 *     6 = vwsub.vx e8 sc_0x00: CSR side-effect
 *     7 = vwsub.vx e8 sc_0x80: result
 *     8 = vwsub.vx e8 sc_0x80: CSR side-effect
 *     9 = vwsub.vx e16 sc_0x0001: result
 *    10 = vwsub.vx e16 sc_0x0001: CSR side-effect
 *    11 = vwsub.vx e16 sc_0xffff: result
 *    12 = vwsub.vx e16 sc_0xffff: CSR side-effect
 *    13 = vwsub.vx e16 sc_0x0000: result
 *    14 = vwsub.vx e16 sc_0x0000: CSR side-effect
 *    15 = vwsub.vx e16 sc_0x8000: result
 *    16 = vwsub.vx e16 sc_0x8000: CSR side-effect
 *    17 = vwsub.vx e32 sc_0x00000001: result
 *    18 = vwsub.vx e32 sc_0x00000001: CSR side-effect
 *    19 = vwsub.vx e32 sc_0xffffffff: result
 *    20 = vwsub.vx e32 sc_0xffffffff: CSR side-effect
 *    21 = vwsub.vx e32 sc_0x00000000: result
 *    22 = vwsub.vx e32 sc_0x00000000: CSR side-effect
 *    23 = vwsub.vx e32 sc_0x80000000: result
 *    24 = vwsub.vx e32 sc_0x80000000: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vwsub.vx SEW=8 sc_0x01 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    li a0, 0x01
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8

    /* Test 3-4: vwsub.vx SEW=8 sc_0xff */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    li a0, 0xff
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 3
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8

    /* Test 5-6: vwsub.vx SEW=8 sc_0x00 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_s2
    vle8.v v16, (t1)
    li a0, 0x00
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 5
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8

    /* Test 7-8: vwsub.vx SEW=8 sc_0x80 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_s2
    vle8.v v16, (t1)
    li a0, 0x80
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 7
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8

    /* Test 9-10: vwsub.vx SEW=16 sc_0x0001 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_s2
    vle16.v v16, (t1)
    li a0, 0x0001
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 9
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16

    /* Test 11-12: vwsub.vx SEW=16 sc_0xffff */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_s2
    vle16.v v16, (t1)
    li a0, 0xffff
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 11
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16

    /* Test 13-14: vwsub.vx SEW=16 sc_0x0000 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc13_s2
    vle16.v v16, (t1)
    li a0, 0x0000
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 14
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 13
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16

    /* Test 15-16: vwsub.vx SEW=16 sc_0x8000 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc15_s2
    vle16.v v16, (t1)
    li a0, 0x8000
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 16
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 15
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 16

    /* Test 17-18: vwsub.vx SEW=32 sc_0x00000001 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc17_s2
    vle32.v v16, (t1)
    li a0, 0x00000001
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 17
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 32

    /* Test 19-20: vwsub.vx SEW=32 sc_0xffffffff */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc19_s2
    vle32.v v16, (t1)
    li a0, 0xffffffff
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 20
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 19
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 32

    /* Test 21-22: vwsub.vx SEW=32 sc_0x00000000 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc21_s2
    vle32.v v16, (t1)
    li a0, 0x00000000
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 22
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 21
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc21_exp, 32

    /* Test 23-24: vwsub.vx SEW=32 sc_0x80000000 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc23_s2
    vle32.v v16, (t1)
    li a0, 0x80000000
    SAVE_CSRS
    vwsub.vx v8, v16, a0
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 23
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc23_exp, 32

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0xff, 0x7f
tc1_exp:
    .half 0x0000, 0x0001, 0xfffe, 0x007e
.align 1
tc3_s2:
    .byte 0x01, 0x02, 0xff, 0x7f
tc3_exp:
    .half 0x0002, 0x0003, 0x0000, 0x0080
.align 1
tc5_s2:
    .byte 0x01, 0x02, 0xff, 0x7f
tc5_exp:
    .half 0x0001, 0x0002, 0xffff, 0x007f
.align 1
tc7_s2:
    .byte 0x01, 0x02, 0xff, 0x7f
tc7_exp:
    .half 0x0081, 0x0082, 0x007f, 0x00ff
.align 2
tc9_s2:
    .half 0x0001, 0x0002, 0xffff, 0x7fff
tc9_exp:
    .word 0x00000000, 0x00000001, 0xfffffffe, 0x00007ffe
.align 2
tc11_s2:
    .half 0x0001, 0x0002, 0xffff, 0x7fff
tc11_exp:
    .word 0x00000002, 0x00000003, 0x00000000, 0x00008000
.align 2
tc13_s2:
    .half 0x0001, 0x0002, 0xffff, 0x7fff
tc13_exp:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x00007fff
.align 2
tc15_s2:
    .half 0x0001, 0x0002, 0xffff, 0x7fff
tc15_exp:
    .word 0x00008001, 0x00008002, 0x00007fff, 0x0000ffff
.align 3
tc17_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x7fffffff
tc17_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0xfffffffffffffffe, 0x000000007ffffffe
.align 3
tc19_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x7fffffff
tc19_exp:
    .dword 0x0000000000000002, 0x0000000000000003, 0x0000000000000000, 0x0000000080000000
.align 3
tc21_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x7fffffff
tc21_exp:
    .dword 0x0000000000000001, 0x0000000000000002, 0xffffffffffffffff, 0x000000007fffffff
.align 3
tc23_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x7fffffff
tc23_exp:
    .dword 0x0000000080000001, 0x0000000080000002, 0x000000007fffffff, 0x00000000ffffffff

.align 4
result_buf:  .space 256
witness_buf: .space 256

