/* Auto-generated test for vlse32.v
 * Strided load of 32-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlse32.v stride=8: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_data
    li t2, 8
    SAVE_CSRS
    vlse32.v v8, (t1), t2
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_data:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004, 0x00000005, 0x00000006, 0x00000007, 0x00000008
tc1_exp:
    .word 0x00000001, 0x00000003, 0x00000005, 0x00000007

.align 4
result_buf:  .space 256
witness_buf: .space 256

