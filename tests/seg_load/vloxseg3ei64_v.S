/* Auto-generated test for vloxseg3ei64.v
 * Indexed segment load nf=3 eew=64
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vloxseg3ei64.v: field 0 result
 *     2 = vloxseg3ei64.v: field 1 result
 *     3 = vloxseg3ei64.v: field 2 result
 *     4 = vloxseg3ei64.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_idx
    vle64.v v16, (t1)
    la t1, tc1_src
    SAVE_CSRS
    vloxseg3ei64.v v8, (t1), v16
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
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_src:
    .dword 0x0000000000000001, 0x0000000000000005, 0x0000000000000009, 0x0000000000000002, 0x0000000000000006, 0x000000000000000a, 0x0000000000000003, 0x0000000000000007, 0x000000000000000b, 0x0000000000000004, 0x0000000000000008, 0x000000000000000c
tc1_idx:
    .dword 0x0000000000000000, 0x0000000000000018, 0x0000000000000030, 0x0000000000000048
tc1_exp_f0:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc1_exp_f1:
    .dword 0x0000000000000005, 0x0000000000000006, 0x0000000000000007, 0x0000000000000008
tc1_exp_f2:
    .dword 0x0000000000000009, 0x000000000000000a, 0x000000000000000b, 0x000000000000000c

.align 4
result_buf:  .space 256
witness_buf: .space 256

