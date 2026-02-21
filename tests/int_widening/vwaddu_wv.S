/* Auto-generated test for vwaddu.wv
 * Widening-W vwaddu.wv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwaddu.wv e8 basic: result
 *     2 = vwaddu.wv e8 basic: CSR side-effect
 *     3 = vwaddu.wv e8 max: result
 *     4 = vwaddu.wv e8 max: CSR side-effect
 *     5 = vwaddu.wv e8 mixed: result
 *     6 = vwaddu.wv e8 mixed: CSR side-effect
 *     7 = vwaddu.wv e16 basic: result
 *     8 = vwaddu.wv e16 basic: CSR side-effect
 *     9 = vwaddu.wv e16 max: result
 *    10 = vwaddu.wv e16 max: CSR side-effect
 *    11 = vwaddu.wv e16 mixed: result
 *    12 = vwaddu.wv e16 mixed: CSR side-effect
 *    13 = vwaddu.wv e32 basic: result
 *    14 = vwaddu.wv e32 basic: CSR side-effect
 *    15 = vwaddu.wv e32 max: result
 *    16 = vwaddu.wv e32 max: CSR side-effect
 *    17 = vwaddu.wv e32 mixed: result
 *    18 = vwaddu.wv e32 mixed: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vwaddu.wv SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc1_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwaddu.wv v8, v16, v20
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8

    /* Test 3-4: vwaddu.wv SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc3_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwaddu.wv v8, v16, v20
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 3
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8

    /* Test 5-6: vwaddu.wv SEW=8 mixed */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc5_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwaddu.wv v8, v16, v20
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 5
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8

    /* Test 7-8: vwaddu.wv SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc7_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwaddu.wv v8, v16, v20
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 7
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 16

    /* Test 9-10: vwaddu.wv SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc9_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwaddu.wv v8, v16, v20
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 9
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16

    /* Test 11-12: vwaddu.wv SEW=16 mixed */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc11_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwaddu.wv v8, v16, v20
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 11
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16

    /* Test 13-14: vwaddu.wv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc13_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwaddu.wv v8, v16, v20
    SET_TEST_NUM 14
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 13
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 32

    /* Test 15-16: vwaddu.wv SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc15_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc15_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwaddu.wv v8, v16, v20
    SET_TEST_NUM 16
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 15
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 32

    /* Test 17-18: vwaddu.wv SEW=32 mixed */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc17_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc17_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwaddu.wv v8, v16, v20
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 17
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 32

    PASS_TEST

.data
.align 1
tc1_s2w:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc1_s1:
    .byte 0x05, 0x06, 0x07, 0x08
tc1_exp:
    .half 0x0006, 0x0008, 0x000a, 0x000c
.align 1
tc3_s2w:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc3_s1:
    .byte 0xff, 0xff, 0xff, 0xff
tc3_exp:
    .half 0x00fe, 0x00fe, 0x00fe, 0x00fe
.align 1
tc5_s2w:
    .half 0x0000, 0x0001, 0x7fff, 0xffff
tc5_s1:
    .byte 0x00, 0x01, 0xff, 0xff
tc5_exp:
    .half 0x0000, 0x0002, 0x80fe, 0x00fe
.align 2
tc7_s2w:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc7_s1:
    .half 0x0005, 0x0006, 0x0007, 0x0008
tc7_exp:
    .word 0x00000006, 0x00000008, 0x0000000a, 0x0000000c
.align 2
tc9_s2w:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc9_s1:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc9_exp:
    .word 0x0000fffe, 0x0000fffe, 0x0000fffe, 0x0000fffe
.align 2
tc11_s2w:
    .word 0x00000000, 0x00000001, 0x7fffffff, 0xffffffff
tc11_s1:
    .half 0x0000, 0x0001, 0xffff, 0xffff
tc11_exp:
    .word 0x00000000, 0x00000002, 0x8000fffe, 0x0000fffe
.align 3
tc13_s2w:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc13_s1:
    .word 0x00000005, 0x00000006, 0x00000007, 0x00000008
tc13_exp:
    .dword 0x0000000000000006, 0x0000000000000008, 0x000000000000000a, 0x000000000000000c
.align 3
tc15_s2w:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
tc15_s1:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc15_exp:
    .dword 0x00000000fffffffe, 0x00000000fffffffe, 0x00000000fffffffe, 0x00000000fffffffe
.align 3
tc17_s2w:
    .dword 0x0000000000000000, 0x0000000000000001, 0x7fffffffffffffff, 0xffffffffffffffff
tc17_s1:
    .word 0x00000000, 0x00000001, 0xffffffff, 0xffffffff
tc17_exp:
    .dword 0x0000000000000000, 0x0000000000000002, 0x80000000fffffffe, 0x00000000fffffffe

.align 4
result_buf:  .space 256
witness_buf: .space 256

