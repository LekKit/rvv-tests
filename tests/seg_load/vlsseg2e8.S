/* Auto-generated test for vlsseg2e8.v
 * Strided segment load nf=2 eew=8
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlsseg2e8.v: field 0 result
 *     2 = vlsseg2e8.v: field 1 result
 *     3 = vlsseg2e8.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    li t2, 2
    SAVE_CSRS
    vlsseg2e8.v v8, (t1), t2
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp_f0, 4
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v9, (t1)
    CHECK_MEM result_buf, tc1_exp_f1, 4
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_src:
    .byte 0x01, 0x05, 0x02, 0x06, 0x03, 0x07, 0x04, 0x08
tc1_exp_f0:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_exp_f1:
    .byte 0x05, 0x06, 0x07, 0x08

.align 4
result_buf:  .space 256
witness_buf: .space 256

