/* Auto-generated test for vlsseg2e32.v
 * Strided segment load nf=2 eew=32
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlsseg2e32.v: field 0 result
 *     2 = vlsseg2e32.v: field 1 result
 *     3 = vlsseg2e32.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    li t2, 8
    SAVE_CSRS
    vlsseg2e32.v v8, (t1), t2
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp_f0, 16
    SET_TEST_NUM 2
    la t1, result_buf
    vse32.v v9, (t1)
    CHECK_MEM result_buf, tc1_exp_f1, 16
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0x00000001, 0x00000005, 0x00000002, 0x00000006, 0x00000003, 0x00000007, 0x00000004, 0x00000008
tc1_exp_f0:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc1_exp_f1:
    .word 0x00000005, 0x00000006, 0x00000007, 0x00000008

.align 4
result_buf:  .space 256
witness_buf: .space 256

