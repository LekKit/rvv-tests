/* Auto-generated test for vsoxseg5ei16.v
 * Indexed segment store nf=5 eew=16
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsoxseg5ei16.v: memory result
 *     2 = vsoxseg5ei16.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_f0
    vle16.v v8, (t1)
    la t1, tc1_f1
    vle16.v v9, (t1)
    la t1, tc1_f2
    vle16.v v10, (t1)
    la t1, tc1_f3
    vle16.v v11, (t1)
    la t1, tc1_f4
    vle16.v v12, (t1)
    la t1, tc1_idx
    vle16.v v16, (t1)
    SAVE_CSRS
    la t1, result_buf
    vsoxseg5ei16.v v8, (t1), v16
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 40
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_f0:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc1_f1:
    .half 0x0005, 0x0006, 0x0007, 0x0008
tc1_f2:
    .half 0x0009, 0x000a, 0x000b, 0x000c
tc1_f3:
    .half 0x000d, 0x000e, 0x000f, 0x0010
tc1_f4:
    .half 0x0011, 0x0012, 0x0013, 0x0014
tc1_idx:
    .half 0x0000, 0x000a, 0x0014, 0x001e
tc1_exp:
    .half 0x0001, 0x0005, 0x0009, 0x000d, 0x0011, 0x0002, 0x0006, 0x000a, 0x000e, 0x0012, 0x0003, 0x0007, 0x000b, 0x000f, 0x0013, 0x0004, 0x0008, 0x000c, 0x0010, 0x0014

.align 4
result_buf:  .space 256
witness_buf: .space 256

