/* Auto-generated test for vwsub.wx
 * Widening-W vwsub.wx
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwsub.wx e8 sc_0x01: result
 *     2 = vwsub.wx e8 sc_0xff: result
 *     3 = vwsub.wx e8 sc_0x00: result
 *     4 = vwsub.wx e16 sc_0x0001: result
 *     5 = vwsub.wx e16 sc_0xffff: result
 *     6 = vwsub.wx e16 sc_0x0000: result
 *     7 = vwsub.wx e32 sc_0x00000001: result
 *     8 = vwsub.wx e32 sc_0xffffffff: result
 *     9 = vwsub.wx e32 sc_0x00000000: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc1_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    li a0, 0x01
    SAVE_CSRS
    vwsub.wx v8, v16, a0
    SET_TEST_NUM 1
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8
    CHECK_VSTART_ZERO

    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc2_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    li a0, 0xff
    SAVE_CSRS
    vwsub.wx v8, v16, a0
    SET_TEST_NUM 2
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_VSTART_ZERO

    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc3_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    li a0, 0x00
    SAVE_CSRS
    vwsub.wx v8, v16, a0
    SET_TEST_NUM 3
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    CHECK_VSTART_ZERO

    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc4_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    li a0, 0x0001
    SAVE_CSRS
    vwsub.wx v8, v16, a0
    SET_TEST_NUM 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 16
    CHECK_VSTART_ZERO

    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc5_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    li a0, 0xffff
    SAVE_CSRS
    vwsub.wx v8, v16, a0
    SET_TEST_NUM 5
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 16
    CHECK_VSTART_ZERO

    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc6_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    li a0, 0x0000
    SAVE_CSRS
    vwsub.wx v8, v16, a0
    SET_TEST_NUM 6
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 16
    CHECK_VSTART_ZERO

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc7_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    li a0, 0x00000001
    SAVE_CSRS
    vwsub.wx v8, v16, a0
    SET_TEST_NUM 7
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 32
    CHECK_VSTART_ZERO

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc8_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    li a0, 0xffffffff
    SAVE_CSRS
    vwsub.wx v8, v16, a0
    SET_TEST_NUM 8
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 32
    CHECK_VSTART_ZERO

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc9_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    li a0, 0x00000000
    SAVE_CSRS
    vwsub.wx v8, v16, a0
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
    .half 0x0001, 0x0002, 0x7fff, 0xffff
tc1_exp:
    .half 0x0000, 0x0001, 0x7ffe, 0xfffe
.align 1
tc2_s2w:
    .half 0x0001, 0x0002, 0x7fff, 0xffff
tc2_exp:
    .half 0x0002, 0x0003, 0x8000, 0x0000
.align 1
tc3_s2w:
    .half 0x0001, 0x0002, 0x7fff, 0xffff
tc3_exp:
    .half 0x0001, 0x0002, 0x7fff, 0xffff
.align 2
tc4_s2w:
    .word 0x00000001, 0x00000002, 0x7fffffff, 0xffffffff
tc4_exp:
    .word 0x00000000, 0x00000001, 0x7ffffffe, 0xfffffffe
.align 2
tc5_s2w:
    .word 0x00000001, 0x00000002, 0x7fffffff, 0xffffffff
tc5_exp:
    .word 0x00000002, 0x00000003, 0x80000000, 0x00000000
.align 2
tc6_s2w:
    .word 0x00000001, 0x00000002, 0x7fffffff, 0xffffffff
tc6_exp:
    .word 0x00000001, 0x00000002, 0x7fffffff, 0xffffffff
.align 3
tc7_s2w:
    .dword 0x0000000000000001, 0x0000000000000002, 0x7fffffffffffffff, 0xffffffffffffffff
tc7_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0x7ffffffffffffffe, 0xfffffffffffffffe
.align 3
tc8_s2w:
    .dword 0x0000000000000001, 0x0000000000000002, 0x7fffffffffffffff, 0xffffffffffffffff
tc8_exp:
    .dword 0x0000000000000002, 0x0000000000000003, 0x8000000000000000, 0x0000000000000000
.align 3
tc9_s2w:
    .dword 0x0000000000000001, 0x0000000000000002, 0x7fffffffffffffff, 0xffffffffffffffff
tc9_exp:
    .dword 0x0000000000000001, 0x0000000000000002, 0x7fffffffffffffff, 0xffffffffffffffff

.align 4
result_buf:  .space 256
witness_buf: .space 256

