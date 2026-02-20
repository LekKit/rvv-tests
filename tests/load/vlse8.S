/* Auto-generated test for vlse8.v
 * Strided load of 8-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlse8.v stride=2: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_data
    li t2, 2
    SAVE_CSRS
    vlse8.v v8, (t1), t2
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .byte 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08
tc1_exp:
    .byte 0x01, 0x03, 0x05, 0x07

.align 4
result_buf:  .space 256
witness_buf: .space 256

