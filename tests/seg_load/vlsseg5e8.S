/* Auto-generated test for vlsseg5e8.v
 * Strided segment load nf=5 eew=8
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlsseg5e8.v: field 0 result
 *     2 = vlsseg5e8.v: field 1 result
 *     3 = vlsseg5e8.v: field 2 result
 *     4 = vlsseg5e8.v: field 3 result
 *     5 = vlsseg5e8.v: field 4 result
 *     6 = vlsseg5e8.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    li t2, 5
    SAVE_CSRS
    vlsseg5e8.v v8, (t1), t2
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
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_src:
    .byte 0x01, 0x05, 0x09, 0x0d, 0x11, 0x02, 0x06, 0x0a, 0x0e, 0x12, 0x03, 0x07, 0x0b, 0x0f, 0x13, 0x04, 0x08, 0x0c, 0x10, 0x14
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

.align 4
result_buf:  .space 256
witness_buf: .space 256

