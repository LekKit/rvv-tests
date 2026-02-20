/* Auto-generated test for vwadd.wv
 * Widening-W vwadd.wv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwadd.wv e8 basic: result
 *     2 = vwadd.wv e8 max: result
 *     3 = vwadd.wv e8 mixed: result
 *     4 = vwadd.wv e16 basic: result
 *     5 = vwadd.wv e16 max: result
 *     6 = vwadd.wv e16 mixed: result
 *     7 = vwadd.wv e32 basic: result
 *     8 = vwadd.wv e32 max: result
 *     9 = vwadd.wv e32 mixed: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vwadd.wv SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc1_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwadd.wv v8, v16, v20
    SET_TEST_NUM 1
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8
    CHECK_VSTART_ZERO

    /* Test 2: vwadd.wv SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc2_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwadd.wv v8, v16, v20
    SET_TEST_NUM 2
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_VSTART_ZERO

    /* Test 3: vwadd.wv SEW=8 mixed */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc3_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwadd.wv v8, v16, v20
    SET_TEST_NUM 3
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    CHECK_VSTART_ZERO

    /* Test 4: vwadd.wv SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc4_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc4_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwadd.wv v8, v16, v20
    SET_TEST_NUM 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 16
    CHECK_VSTART_ZERO

    /* Test 5: vwadd.wv SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc5_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwadd.wv v8, v16, v20
    SET_TEST_NUM 5
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 16
    CHECK_VSTART_ZERO

    /* Test 6: vwadd.wv SEW=16 mixed */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc6_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwadd.wv v8, v16, v20
    SET_TEST_NUM 6
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 16
    CHECK_VSTART_ZERO

    /* Test 7: vwadd.wv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc7_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc7_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwadd.wv v8, v16, v20
    SET_TEST_NUM 7
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 32
    CHECK_VSTART_ZERO

    /* Test 8: vwadd.wv SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc8_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc8_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwadd.wv v8, v16, v20
    SET_TEST_NUM 8
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 32
    CHECK_VSTART_ZERO

    /* Test 9: vwadd.wv SEW=32 mixed */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc9_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwadd.wv v8, v16, v20
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
    .half 0x0006, 0x0008, 0x000a, 0x000c
.align 1
tc2_s2w:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc2_s1:
    .byte 0xff, 0xff, 0xff, 0xff
tc2_exp:
    .half 0xfffe, 0xfffe, 0xfffe, 0xfffe
.align 1
tc3_s2w:
    .half 0x0000, 0x0001, 0x7fff, 0xffff
tc3_s1:
    .byte 0x00, 0x01, 0xff, 0xff
tc3_exp:
    .half 0x0000, 0x0002, 0x7ffe, 0xfffe
.align 2
tc4_s2w:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc4_s1:
    .half 0x0005, 0x0006, 0x0007, 0x0008
tc4_exp:
    .word 0x00000006, 0x00000008, 0x0000000a, 0x0000000c
.align 2
tc5_s2w:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc5_s1:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc5_exp:
    .word 0xfffffffe, 0xfffffffe, 0xfffffffe, 0xfffffffe
.align 2
tc6_s2w:
    .word 0x00000000, 0x00000001, 0x7fffffff, 0xffffffff
tc6_s1:
    .half 0x0000, 0x0001, 0xffff, 0xffff
tc6_exp:
    .word 0x00000000, 0x00000002, 0x7ffffffe, 0xfffffffe
.align 3
tc7_s2w:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc7_s1:
    .word 0x00000005, 0x00000006, 0x00000007, 0x00000008
tc7_exp:
    .dword 0x0000000000000006, 0x0000000000000008, 0x000000000000000a, 0x000000000000000c
.align 3
tc8_s2w:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
tc8_s1:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc8_exp:
    .dword 0xfffffffffffffffe, 0xfffffffffffffffe, 0xfffffffffffffffe, 0xfffffffffffffffe
.align 3
tc9_s2w:
    .dword 0x0000000000000000, 0x0000000000000001, 0x7fffffffffffffff, 0xffffffffffffffff
tc9_s1:
    .word 0x00000000, 0x00000001, 0xffffffff, 0xffffffff
tc9_exp:
    .dword 0x0000000000000000, 0x0000000000000002, 0x7ffffffffffffffe, 0xfffffffffffffffe

.align 4
result_buf:  .space 256
witness_buf: .space 256

