/* Auto-generated test for vlseg6e32.v
 * Unit-stride segment load nf=6 eew=32
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlseg6e32.v: field 0 result
 *     2 = vlseg6e32.v: field 1 result
 *     3 = vlseg6e32.v: field 2 result
 *     4 = vlseg6e32.v: field 3 result
 *     5 = vlseg6e32.v: field 4 result
 *     6 = vlseg6e32.v: field 5 result
 *     7 = vlseg6e32.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    SAVE_CSRS
    vlseg6e32.v v8, (t1)
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp_f0, 16
    SET_TEST_NUM 2
    la t1, result_buf
    vse32.v v9, (t1)
    CHECK_MEM result_buf, tc1_exp_f1, 16
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v10, (t1)
    CHECK_MEM result_buf, tc1_exp_f2, 16
    SET_TEST_NUM 4
    la t1, result_buf
    vse32.v v11, (t1)
    CHECK_MEM result_buf, tc1_exp_f3, 16
    SET_TEST_NUM 5
    la t1, result_buf
    vse32.v v12, (t1)
    CHECK_MEM result_buf, tc1_exp_f4, 16
    SET_TEST_NUM 6
    la t1, result_buf
    vse32.v v13, (t1)
    CHECK_MEM result_buf, tc1_exp_f5, 16
    SET_TEST_NUM 7
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0x00000001, 0x00000005, 0x00000009, 0x0000000d, 0x00000011, 0x00000015, 0x00000002, 0x00000006, 0x0000000a, 0x0000000e, 0x00000012, 0x00000016, 0x00000003, 0x00000007, 0x0000000b, 0x0000000f, 0x00000013, 0x00000017, 0x00000004, 0x00000008, 0x0000000c, 0x00000010, 0x00000014, 0x00000018
tc1_exp_f0:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc1_exp_f1:
    .word 0x00000005, 0x00000006, 0x00000007, 0x00000008
tc1_exp_f2:
    .word 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c
tc1_exp_f3:
    .word 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010
tc1_exp_f4:
    .word 0x00000011, 0x00000012, 0x00000013, 0x00000014
tc1_exp_f5:
    .word 0x00000015, 0x00000016, 0x00000017, 0x00000018

.align 4
result_buf:  .space 256
witness_buf: .space 256

