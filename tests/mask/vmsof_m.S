/* Auto-generated test for vmsof.m
 * Set only first set mask bit
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmsof.m none: result
 *     2 = vmsof.m bit0: result
 *     3 = vmsof.m bit1: result
 *     4 = vmsof.m bit2: result
 *     5 = vmsof.m bit3: result
 *     6 = vmsof.m multi: result
 *     7 = vmsof.m all: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vmsof.m v8, v16
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
    la t1, tc2_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vmsof.m v8, v16
    SET_TEST_NUM 2
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 1
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vmsof.m v8, v16
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 2
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vmsof.m v8, v16
    SET_TEST_NUM 4
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 4
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vmsof.m v8, v16
    SET_TEST_NUM 5
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 8
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc6_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vmsof.m v8, v16
    SET_TEST_NUM 6
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 2
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vmsof.m v8, v16
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 1
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_mask:
    .byte 0
.align 1
tc2_mask:
    .byte 1
.align 1
tc3_mask:
    .byte 2
.align 1
tc4_mask:
    .byte 4
.align 1
tc5_mask:
    .byte 8
.align 1
tc6_mask:
    .byte 10
.align 1
tc7_mask:
    .byte 15

.align 4
result_buf:  .space 256
witness_buf: .space 256

