/* Auto-generated test for vfirst.m
 * Find first set mask bit
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfirst.m none: wrong index
 *     2 = vfirst.m bit0: wrong index
 *     3 = vfirst.m bit1: wrong index
 *     4 = vfirst.m bit2: wrong index
 *     5 = vfirst.m bit3: wrong index
 *     6 = vfirst.m multi: wrong index
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vfirst.m t2, v16
    SET_TEST_NUM 1
    li t3, -1
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vfirst.m t2, v16
    SET_TEST_NUM 2
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vfirst.m t2, v16
    SET_TEST_NUM 3
    li t3, 1
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vfirst.m t2, v16
    SET_TEST_NUM 4
    li t3, 2
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vfirst.m t2, v16
    SET_TEST_NUM 5
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc6_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vfirst.m t2, v16
    SET_TEST_NUM 6
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

.align 4
result_buf:  .space 256
witness_buf: .space 256

