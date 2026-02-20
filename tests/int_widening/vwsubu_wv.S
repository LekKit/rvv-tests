/* Auto-generated test for vwsubu.wv
 * Widening-W vwsubu.wv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwsubu.wv e8 basic: result
 *     2 = vwsubu.wv e8 max: result
 *     3 = vwsubu.wv e8 mixed: result
 *     4 = vwsubu.wv e16 basic: result
 *     5 = vwsubu.wv e16 max: result
 *     6 = vwsubu.wv e16 mixed: result
 *     7 = vwsubu.wv e32 basic: result
 *     8 = vwsubu.wv e32 max: result
 *     9 = vwsubu.wv e32 mixed: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vwsubu.wv SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc1_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwsubu.wv v8, v16, v20
    SET_TEST_NUM 1
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8
    CHECK_VSTART_ZERO

    /* Test 2: vwsubu.wv SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc2_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwsubu.wv v8, v16, v20
    SET_TEST_NUM 2
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_VSTART_ZERO

    /* Test 3: vwsubu.wv SEW=8 mixed */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc3_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwsubu.wv v8, v16, v20
    SET_TEST_NUM 3
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    CHECK_VSTART_ZERO

    /* Test 4: vwsubu.wv SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc4_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc4_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwsubu.wv v8, v16, v20
    SET_TEST_NUM 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 16
    CHECK_VSTART_ZERO

    /* Test 5: vwsubu.wv SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc5_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwsubu.wv v8, v16, v20
    SET_TEST_NUM 5
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 16
    CHECK_VSTART_ZERO

    /* Test 6: vwsubu.wv SEW=16 mixed */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc6_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwsubu.wv v8, v16, v20
    SET_TEST_NUM 6
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 16
    CHECK_VSTART_ZERO

    /* Test 7: vwsubu.wv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc7_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc7_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwsubu.wv v8, v16, v20
    SET_TEST_NUM 7
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 32
    CHECK_VSTART_ZERO

    /* Test 8: vwsubu.wv SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc8_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc8_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwsubu.wv v8, v16, v20
    SET_TEST_NUM 8
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 32
    CHECK_VSTART_ZERO

    /* Test 9: vwsubu.wv SEW=32 mixed */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc9_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwsubu.wv v8, v16, v20
    SET_TEST_NUM 9
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 32
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 1
tc1_s2w:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc1_s1:
    .byte 0x05, 0x06, 0x07, 0x08
tc1_exp:
    .half 0xfffc, 0xfffc, 0xfffc, 0xfffc
.align 1
tc2_s2w:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc2_s1:
    .byte 0xff, 0xff, 0xff, 0xff
tc2_exp:
    .half 0xff00, 0xff00, 0xff00, 0xff00
.align 1
tc3_s2w:
    .half 0x0000, 0x0001, 0x7fff, 0xffff
tc3_s1:
    .byte 0x00, 0x01, 0xff, 0xff
tc3_exp:
    .half 0x0000, 0x0000, 0x7f00, 0xff00
.align 2
tc4_s2w:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc4_s1:
    .half 0x0005, 0x0006, 0x0007, 0x0008
tc4_exp:
    .word 0xfffffffc, 0xfffffffc, 0xfffffffc, 0xfffffffc
.align 2
tc5_s2w:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc5_s1:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc5_exp:
    .word 0xffff0000, 0xffff0000, 0xffff0000, 0xffff0000
.align 2
tc6_s2w:
    .word 0x00000000, 0x00000001, 0x7fffffff, 0xffffffff
tc6_s1:
    .half 0x0000, 0x0001, 0xffff, 0xffff
tc6_exp:
    .word 0x00000000, 0x00000000, 0x7fff0000, 0xffff0000
.align 3
tc7_s2w:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc7_s1:
    .word 0x00000005, 0x00000006, 0x00000007, 0x00000008
tc7_exp:
    .dword 0xfffffffffffffffc, 0xfffffffffffffffc, 0xfffffffffffffffc, 0xfffffffffffffffc
.align 3
tc8_s2w:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
tc8_s1:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc8_exp:
    .dword 0xffffffff00000000, 0xffffffff00000000, 0xffffffff00000000, 0xffffffff00000000
.align 3
tc9_s2w:
    .dword 0x0000000000000000, 0x0000000000000001, 0x7fffffffffffffff, 0xffffffffffffffff
tc9_s1:
    .word 0x00000000, 0x00000001, 0xffffffff, 0xffffffff
tc9_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x7fffffff00000000, 0xffffffff00000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

