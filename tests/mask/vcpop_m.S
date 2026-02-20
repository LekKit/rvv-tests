/* Auto-generated test for vcpop.m
 * Mask population count
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vcpop.m none: wrong count
 *     2 = vcpop.m one: wrong count
 *     3 = vcpop.m two: wrong count
 *     4 = vcpop.m all: wrong count
 *     5 = vcpop.m three: wrong count
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vcpop.m none */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vcpop.m t2, v16
    SET_TEST_NUM 1
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    /* Test 2: vcpop.m one */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vcpop.m t2, v16
    SET_TEST_NUM 2
    li t3, 1
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    /* Test 3: vcpop.m two */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vcpop.m t2, v16
    SET_TEST_NUM 3
    li t3, 2
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    /* Test 4: vcpop.m all */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vcpop.m t2, v16
    SET_TEST_NUM 4
    li t3, 4
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    /* Test 5: vcpop.m three */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_mask
    vlm.v v16, (t1)
    SAVE_CSRS
    vcpop.m t2, v16
    SET_TEST_NUM 5
    li t3, 3
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
    .byte 5
.align 1
tc4_mask:
    .byte 15
.align 1
tc5_mask:
    .byte 14

.align 4
result_buf:  .space 256
witness_buf: .space 256

