/* Auto-generated test for vlseg3e16ff.v
 * Fault-first segment load nf=3 eew=16
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlseg3e16ff.v: field 0 result
 *     2 = vlseg3e16ff.v: field 1 result
 *     3 = vlseg3e16ff.v: field 2 result
 *     4 = vlseg3e16ff.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_src
    SAVE_CSRS
    vlseg3e16ff.v v8, (t1)
    SET_TEST_NUM 1
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp_f0, 8
    SET_TEST_NUM 2
    la t1, result_buf
    vse16.v v9, (t1)
    CHECK_MEM result_buf, tc1_exp_f1, 8
    SET_TEST_NUM 3
    la t1, result_buf
    vse16.v v10, (t1)
    CHECK_MEM result_buf, tc1_exp_f2, 8
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_src:
    .half 0x0001, 0x0005, 0x0009, 0x0002, 0x0006, 0x000a, 0x0003, 0x0007, 0x000b, 0x0004, 0x0008, 0x000c
tc1_exp_f0:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc1_exp_f1:
    .half 0x0005, 0x0006, 0x0007, 0x0008
tc1_exp_f2:
    .half 0x0009, 0x000a, 0x000b, 0x000c

.align 4
result_buf:  .space 256
witness_buf: .space 256

