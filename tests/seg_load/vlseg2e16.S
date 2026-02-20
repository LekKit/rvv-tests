/* Auto-generated test for vlseg2e16.v
 * Unit-stride segment load nf=2 eew=16
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlseg2e16.v: field 0 result
 *     2 = vlseg2e16.v: field 1 result
 *     3 = vlseg2e16.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_src
    SAVE_CSRS
    vlseg2e16.v v8, (t1)
    SET_TEST_NUM 1
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp_f0, 8
    SET_TEST_NUM 2
    la t1, result_buf
    vse16.v v9, (t1)
    CHECK_MEM result_buf, tc1_exp_f1, 8
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_src:
    .half 0x0001, 0x0005, 0x0002, 0x0006, 0x0003, 0x0007, 0x0004, 0x0008
tc1_exp_f0:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc1_exp_f1:
    .half 0x0005, 0x0006, 0x0007, 0x0008

.align 4
result_buf:  .space 256
witness_buf: .space 256

