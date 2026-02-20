/* Auto-generated test for vloxseg4ei32.v
 * Indexed segment load nf=4 eew=32
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vloxseg4ei32.v: field 0 result
 *     2 = vloxseg4ei32.v: field 1 result
 *     3 = vloxseg4ei32.v: field 2 result
 *     4 = vloxseg4ei32.v: field 3 result
 *     5 = vloxseg4ei32.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_idx
    vle32.v v16, (t1)
    la t1, tc1_src
    SAVE_CSRS
    vloxseg4ei32.v v8, (t1), v16
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
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0x00000001, 0x00000005, 0x00000009, 0x0000000d, 0x00000002, 0x00000006, 0x0000000a, 0x0000000e, 0x00000003, 0x00000007, 0x0000000b, 0x0000000f, 0x00000004, 0x00000008, 0x0000000c, 0x00000010
tc1_idx:
    .word 0x00000000, 0x00000010, 0x00000020, 0x00000030
tc1_exp_f0:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc1_exp_f1:
    .word 0x00000005, 0x00000006, 0x00000007, 0x00000008
tc1_exp_f2:
    .word 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c
tc1_exp_f3:
    .word 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010

.align 4
result_buf:  .space 256
witness_buf: .space 256

