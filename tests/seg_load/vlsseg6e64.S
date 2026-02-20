/* Auto-generated test for vlsseg6e64.v
 * Strided segment load nf=6 eew=64
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlsseg6e64.v: field 0 result
 *     2 = vlsseg6e64.v: field 1 result
 *     3 = vlsseg6e64.v: field 2 result
 *     4 = vlsseg6e64.v: field 3 result
 *     5 = vlsseg6e64.v: field 4 result
 *     6 = vlsseg6e64.v: field 5 result
 *     7 = vlsseg6e64.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_src
    li t2, 48
    SAVE_CSRS
    vlsseg6e64.v v8, (t1), t2
    SET_TEST_NUM 1
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp_f0, 32
    SET_TEST_NUM 2
    la t1, result_buf
    vse64.v v9, (t1)
    CHECK_MEM result_buf, tc1_exp_f1, 32
    SET_TEST_NUM 3
    la t1, result_buf
    vse64.v v10, (t1)
    CHECK_MEM result_buf, tc1_exp_f2, 32
    SET_TEST_NUM 4
    la t1, result_buf
    vse64.v v11, (t1)
    CHECK_MEM result_buf, tc1_exp_f3, 32
    SET_TEST_NUM 5
    la t1, result_buf
    vse64.v v12, (t1)
    CHECK_MEM result_buf, tc1_exp_f4, 32
    SET_TEST_NUM 6
    la t1, result_buf
    vse64.v v13, (t1)
    CHECK_MEM result_buf, tc1_exp_f5, 32
    SET_TEST_NUM 7
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_src:
    .dword 0x0000000000000001, 0x0000000000000005, 0x0000000000000009, 0x000000000000000d, 0x0000000000000011, 0x0000000000000015, 0x0000000000000002, 0x0000000000000006, 0x000000000000000a, 0x000000000000000e, 0x0000000000000012, 0x0000000000000016, 0x0000000000000003, 0x0000000000000007, 0x000000000000000b, 0x000000000000000f, 0x0000000000000013, 0x0000000000000017, 0x0000000000000004, 0x0000000000000008, 0x000000000000000c, 0x0000000000000010, 0x0000000000000014, 0x0000000000000018
tc1_exp_f0:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc1_exp_f1:
    .dword 0x0000000000000005, 0x0000000000000006, 0x0000000000000007, 0x0000000000000008
tc1_exp_f2:
    .dword 0x0000000000000009, 0x000000000000000a, 0x000000000000000b, 0x000000000000000c
tc1_exp_f3:
    .dword 0x000000000000000d, 0x000000000000000e, 0x000000000000000f, 0x0000000000000010
tc1_exp_f4:
    .dword 0x0000000000000011, 0x0000000000000012, 0x0000000000000013, 0x0000000000000014
tc1_exp_f5:
    .dword 0x0000000000000015, 0x0000000000000016, 0x0000000000000017, 0x0000000000000018

.align 4
result_buf:  .space 256
witness_buf: .space 256

