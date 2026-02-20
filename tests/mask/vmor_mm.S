/* Auto-generated test for vmor.mm
 * Mask logical vmor.mm
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmor.mm all_ones: wrong mask result
 *     2 = vmor.mm all_ones: CSR side-effect
 *     3 = vmor.mm all_zeros: wrong mask result
 *     4 = vmor.mm all_zeros: CSR side-effect
 *     5 = vmor.mm mixed_1: wrong mask result
 *     6 = vmor.mm mixed_1: CSR side-effect
 *     7 = vmor.mm mixed_2: wrong mask result
 *     8 = vmor.mm mixed_2: CSR side-effect
 *     9 = vmor.mm one_zero: wrong mask result
 *    10 = vmor.mm one_zero: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vmor.mm all_ones */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_m2
    vlm.v v16, (t1)
    la t1, tc1_m1
    vlm.v v20, (t1)
    SAVE_CSRS
    vmor.mm v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    /* Test 3-4: vmor.mm all_zeros */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_m2
    vlm.v v16, (t1)
    la t1, tc3_m1
    vlm.v v20, (t1)
    SAVE_CSRS
    vmor.mm v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED

    /* Test 5-6: vmor.mm mixed_1 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_m2
    vlm.v v16, (t1)
    la t1, tc5_m1
    vlm.v v20, (t1)
    SAVE_CSRS
    vmor.mm v8, v16, v20
    SET_TEST_NUM 5
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 14
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-8: vmor.mm mixed_2 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_m2
    vlm.v v16, (t1)
    la t1, tc7_m1
    vlm.v v20, (t1)
    SAVE_CSRS
    vmor.mm v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 7
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED

    /* Test 9-10: vmor.mm one_zero */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc9_m2
    vlm.v v16, (t1)
    la t1, tc9_m1
    vlm.v v20, (t1)
    SAVE_CSRS
    vmor.mm v8, v16, v20
    SET_TEST_NUM 9
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_m2:
    .byte 15
tc1_m1:
    .byte 15
.align 1
tc3_m2:
    .byte 0
tc3_m1:
    .byte 0
.align 1
tc5_m2:
    .byte 10
tc5_m1:
    .byte 12
.align 1
tc7_m2:
    .byte 5
tc7_m1:
    .byte 3
.align 1
tc9_m2:
    .byte 15
tc9_m1:
    .byte 0

.align 4
result_buf:  .space 256
witness_buf: .space 256

