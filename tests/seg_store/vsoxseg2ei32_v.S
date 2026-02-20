/* Auto-generated test for vsoxseg2ei32.v
 * Indexed segment store nf=2 eew=32
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsoxseg2ei32.v: memory result
 *     2 = vsoxseg2ei32.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_f0
    vle32.v v8, (t1)
    la t1, tc1_f1
    vle32.v v9, (t1)
    la t1, tc1_idx
    vle32.v v16, (t1)
    SAVE_CSRS
    la t1, result_buf
    vsoxseg2ei32.v v8, (t1), v16
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 32
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_f0:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc1_f1:
    .word 0x00000005, 0x00000006, 0x00000007, 0x00000008
tc1_idx:
    .word 0x00000000, 0x00000008, 0x00000010, 0x00000018
tc1_exp:
    .word 0x00000001, 0x00000005, 0x00000002, 0x00000006, 0x00000003, 0x00000007, 0x00000004, 0x00000008

.align 4
result_buf:  .space 256
witness_buf: .space 256

