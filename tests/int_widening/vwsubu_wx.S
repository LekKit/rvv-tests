/* Auto-generated test for vwsubu.wx
 * Widening-W vwsubu.wx
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwsubu.wx e8 sc_0x01: result
 *     2 = vwsubu.wx e8 sc_0x01: CSR side-effect
 *     3 = vwsubu.wx e8 sc_0xff: result
 *     4 = vwsubu.wx e8 sc_0xff: CSR side-effect
 *     5 = vwsubu.wx e8 sc_0x00: result
 *     6 = vwsubu.wx e8 sc_0x00: CSR side-effect
 *     7 = vwsubu.wx e16 sc_0x0001: result
 *     8 = vwsubu.wx e16 sc_0x0001: CSR side-effect
 *     9 = vwsubu.wx e16 sc_0xffff: result
 *    10 = vwsubu.wx e16 sc_0xffff: CSR side-effect
 *    11 = vwsubu.wx e16 sc_0x0000: result
 *    12 = vwsubu.wx e16 sc_0x0000: CSR side-effect
 *    13 = vwsubu.wx e32 sc_0x00000001: result
 *    14 = vwsubu.wx e32 sc_0x00000001: CSR side-effect
 *    15 = vwsubu.wx e32 sc_0xffffffff: result
 *    16 = vwsubu.wx e32 sc_0xffffffff: CSR side-effect
 *    17 = vwsubu.wx e32 sc_0x00000000: result
 *    18 = vwsubu.wx e32 sc_0x00000000: CSR side-effect
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
    vwsubu.wx v8, v16, a0
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8

    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc3_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    li a0, 0xff
    SAVE_CSRS
    vwsubu.wx v8, v16, a0
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 3
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8

    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc5_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    li a0, 0x00
    SAVE_CSRS
    vwsubu.wx v8, v16, a0
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 5
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8

    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc7_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    li a0, 0x0001
    SAVE_CSRS
    vwsubu.wx v8, v16, a0
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 7
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 16

    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc9_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    li a0, 0xffff
    SAVE_CSRS
    vwsubu.wx v8, v16, a0
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 9
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16

    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc11_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    li a0, 0x0000
    SAVE_CSRS
    vwsubu.wx v8, v16, a0
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 11
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc13_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    li a0, 0x00000001
    SAVE_CSRS
    vwsubu.wx v8, v16, a0
    SET_TEST_NUM 14
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 13
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 32

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc15_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    li a0, 0xffffffff
    SAVE_CSRS
    vwsubu.wx v8, v16, a0
    SET_TEST_NUM 16
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 15
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 32

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc17_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    li a0, 0x00000000
    SAVE_CSRS
    vwsubu.wx v8, v16, a0
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
    .half 0x0001, 0x0002, 0x7fff, 0xffff
tc1_exp:
    .half 0x0000, 0x0001, 0x7ffe, 0xfffe
.align 1
tc3_s2w:
    .half 0x0001, 0x0002, 0x7fff, 0xffff
tc3_exp:
    .half 0xff02, 0xff03, 0x7f00, 0xff00
.align 1
tc5_s2w:
    .half 0x0001, 0x0002, 0x7fff, 0xffff
tc5_exp:
    .half 0x0001, 0x0002, 0x7fff, 0xffff
.align 2
tc7_s2w:
    .word 0x00000001, 0x00000002, 0x7fffffff, 0xffffffff
tc7_exp:
    .word 0x00000000, 0x00000001, 0x7ffffffe, 0xfffffffe
.align 2
tc9_s2w:
    .word 0x00000001, 0x00000002, 0x7fffffff, 0xffffffff
tc9_exp:
    .word 0xffff0002, 0xffff0003, 0x7fff0000, 0xffff0000
.align 2
tc11_s2w:
    .word 0x00000001, 0x00000002, 0x7fffffff, 0xffffffff
tc11_exp:
    .word 0x00000001, 0x00000002, 0x7fffffff, 0xffffffff
.align 3
tc13_s2w:
    .dword 0x0000000000000001, 0x0000000000000002, 0x7fffffffffffffff, 0xffffffffffffffff
tc13_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0x7ffffffffffffffe, 0xfffffffffffffffe
.align 3
tc15_s2w:
    .dword 0x0000000000000001, 0x0000000000000002, 0x7fffffffffffffff, 0xffffffffffffffff
tc15_exp:
    .dword 0xffffffff00000002, 0xffffffff00000003, 0x7fffffff00000000, 0xffffffff00000000
.align 3
tc17_s2w:
    .dword 0x0000000000000001, 0x0000000000000002, 0x7fffffffffffffff, 0xffffffffffffffff
tc17_exp:
    .dword 0x0000000000000001, 0x0000000000000002, 0x7fffffffffffffff, 0xffffffffffffffff

.align 4
result_buf:  .space 256
witness_buf: .space 256

