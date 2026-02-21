/* Auto-generated test for vmnor.mm
 * Mask logical vmnor.mm
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmnor.mm all_ones: wrong mask result
 *     2 = vmnor.mm all_ones: witness changed
 *     3 = vmnor.mm all_ones: CSR side-effect
 *     4 = vmnor.mm all_zeros: wrong mask result
 *     5 = vmnor.mm all_zeros: witness changed
 *     6 = vmnor.mm all_zeros: CSR side-effect
 *     7 = vmnor.mm mixed_1: wrong mask result
 *     8 = vmnor.mm mixed_1: witness changed
 *     9 = vmnor.mm mixed_1: CSR side-effect
 *    10 = vmnor.mm mixed_2: wrong mask result
 *    11 = vmnor.mm mixed_2: witness changed
 *    12 = vmnor.mm mixed_2: CSR side-effect
 *    13 = vmnor.mm one_zero: wrong mask result
 *    14 = vmnor.mm one_zero: witness changed
 *    15 = vmnor.mm one_zero: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vmnor.mm all_ones */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_m2
    vlm.v v16, (t1)
    la t1, tc1_m1
    vlm.v v20, (t1)
    la t1, tc1_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmnor.mm v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 2
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 4
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED

    /* Test 4-6: vmnor.mm all_zeros */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_m2
    vlm.v v16, (t1)
    la t1, tc4_m1
    vlm.v v20, (t1)
    la t1, tc4_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmnor.mm v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 5
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 4
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-9: vmnor.mm mixed_1 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_m2
    vlm.v v16, (t1)
    la t1, tc7_m1
    vlm.v v20, (t1)
    la t1, tc7_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmnor.mm v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 1
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 8
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 4
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED

    /* Test 10-12: vmnor.mm mixed_2 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc10_m2
    vlm.v v16, (t1)
    la t1, tc10_m1
    vlm.v v20, (t1)
    la t1, tc10_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmnor.mm v8, v16, v20
    SET_TEST_NUM 10
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 8
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 11
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 4
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13-15: vmnor.mm one_zero */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc13_m2
    vlm.v v16, (t1)
    la t1, tc13_m1
    vlm.v v20, (t1)
    la t1, tc13_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmnor.mm v8, v16, v20
    SET_TEST_NUM 13
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 14
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 4
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_m2:
    .byte 15
tc1_m1:
    .byte 15
tc1_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc4_m2:
    .byte 0
tc4_m1:
    .byte 0
tc4_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc7_m2:
    .byte 10
tc7_m1:
    .byte 12
tc7_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc10_m2:
    .byte 5
tc10_m1:
    .byte 3
tc10_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc13_m2:
    .byte 15
tc13_m1:
    .byte 0
tc13_w:
    .byte 0xde, 0xad, 0xbe, 0xef

.align 4
result_buf:  .space 256
witness_buf: .space 256

