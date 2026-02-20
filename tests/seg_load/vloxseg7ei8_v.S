/* Auto-generated test for vloxseg7ei8.v
 * Indexed segment load nf=7 eew=8
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vloxseg7ei8.v: field 0 result
 *     2 = vloxseg7ei8.v: field 1 result
 *     3 = vloxseg7ei8.v: field 2 result
 *     4 = vloxseg7ei8.v: field 3 result
 *     5 = vloxseg7ei8.v: field 4 result
 *     6 = vloxseg7ei8.v: field 5 result
 *     7 = vloxseg7ei8.v: field 6 result
 *     8 = vloxseg7ei8.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_idx
    vle8.v v16, (t1)
    la t1, tc1_src
    SAVE_CSRS
    vloxseg7ei8.v v8, (t1), v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp_f0, 4
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v9, (t1)
    CHECK_MEM result_buf, tc1_exp_f1, 4
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v10, (t1)
    CHECK_MEM result_buf, tc1_exp_f2, 4
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v11, (t1)
    CHECK_MEM result_buf, tc1_exp_f3, 4
    SET_TEST_NUM 5
    la t1, result_buf
    vse8.v v12, (t1)
    CHECK_MEM result_buf, tc1_exp_f4, 4
    SET_TEST_NUM 6
    la t1, result_buf
    vse8.v v13, (t1)
    CHECK_MEM result_buf, tc1_exp_f5, 4
    SET_TEST_NUM 7
    la t1, result_buf
    vse8.v v14, (t1)
    CHECK_MEM result_buf, tc1_exp_f6, 4
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_src:
    .byte 0x01, 0x05, 0x09, 0x0d, 0x11, 0x15, 0x19, 0x02, 0x06, 0x0a, 0x0e, 0x12, 0x16, 0x1a, 0x03, 0x07, 0x0b, 0x0f, 0x13, 0x17, 0x1b, 0x04, 0x08, 0x0c, 0x10, 0x14, 0x18, 0x1c
tc1_idx:
    .byte 0x00, 0x07, 0x0e, 0x15
tc1_exp_f0:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_exp_f1:
    .byte 0x05, 0x06, 0x07, 0x08
tc1_exp_f2:
    .byte 0x09, 0x0a, 0x0b, 0x0c
tc1_exp_f3:
    .byte 0x0d, 0x0e, 0x0f, 0x10
tc1_exp_f4:
    .byte 0x11, 0x12, 0x13, 0x14
tc1_exp_f5:
    .byte 0x15, 0x16, 0x17, 0x18
tc1_exp_f6:
    .byte 0x19, 0x1a, 0x1b, 0x1c

.align 4
result_buf:  .space 256
witness_buf: .space 256

