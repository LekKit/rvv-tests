/* Auto-generated test for vssseg2e64.v
 * Strided segment store nf=2 eew=64
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vssseg2e64.v: memory result
 *     2 = vssseg2e64.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_f0
    vle64.v v8, (t1)
    la t1, tc1_f1
    vle64.v v9, (t1)
    SAVE_CSRS
    la t1, result_buf
    li t2, 16
    vssseg2e64.v v8, (t1), t2
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 64
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_f0:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc1_f1:
    .dword 0x0000000000000005, 0x0000000000000006, 0x0000000000000007, 0x0000000000000008
tc1_exp:
    .dword 0x0000000000000001, 0x0000000000000005, 0x0000000000000002, 0x0000000000000006, 0x0000000000000003, 0x0000000000000007, 0x0000000000000004, 0x0000000000000008

.align 4
result_buf:  .space 256
witness_buf: .space 256

