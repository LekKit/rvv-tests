/* Auto-generated test for vmadc.vx
 * Carry-out detection (no carry-in, vx)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmadc.vx e8 basic: result
 *     2 = vmadc.vx e8 overflow: result
 *     3 = vmadc.vx e8 max: result
 *     4 = vmadc.vx e16 basic: result
 *     5 = vmadc.vx e16 overflow: result
 *     6 = vmadc.vx e16 max: result
 *     7 = vmadc.vx e32 basic: result
 *     8 = vmadc.vx e32 overflow: result
 *     9 = vmadc.vx e32 max: result
 *    10 = vmadc.vx e64 basic: result
 *    11 = vmadc.vx e64 overflow: result
 *    12 = vmadc.vx e64 max: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    li a0, 0x0a
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_s2
    vle8.v v16, (t1)
    li a0, 0x01
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 2
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    li a0, 0xff
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc4_s2
    vle16.v v16, (t1)
    li a0, 0x000a
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 4
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_s2
    vle16.v v16, (t1)
    li a0, 0x0001
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 5
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_s2
    vle16.v v16, (t1)
    li a0, 0xffff
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 6
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc7_s2
    vle32.v v16, (t1)
    li a0, 0x0000000a
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc8_s2
    vle32.v v16, (t1)
    li a0, 0x00000001
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 8
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_s2
    vle32.v v16, (t1)
    li a0, 0xffffffff
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 9
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc10_s2
    vle64.v v16, (t1)
    li a0, 0x000000000000000a
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 10
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc11_s2
    vle64.v v16, (t1)
    li a0, 0x0000000000000001
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 11
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc12_s2
    vle64.v v16, (t1)
    li a0, 0xffffffffffffffff
    SAVE_CSRS
    vmadc.vx v8, v16, a0
    SET_TEST_NUM 12
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc2_s2:
    .byte 0xff, 0xff, 0xfe, 0x00
.align 1
tc3_s2:
    .byte 0xff, 0xff, 0xff, 0xff
.align 1
tc4_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc5_s2:
    .half 0xffff, 0xffff, 0xfffe, 0x0000
.align 1
tc6_s2:
    .half 0xffff, 0xffff, 0xffff, 0xffff
.align 2
tc7_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc8_s2:
    .word 0xffffffff, 0xffffffff, 0xfffffffe, 0x00000000
.align 2
tc9_s2:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
.align 3
tc10_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc11_s2:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000
.align 3
tc12_s2:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff

.align 4
result_buf:  .space 256
witness_buf: .space 256

