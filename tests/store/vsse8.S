/* Auto-generated test for vsse8.v
 * Strided store of 8-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsse8.v stride=2: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_data
    vle8.v v8, (t1)
    la t1, result_buf
    li t2, 2
    SAVE_CSRS
    vsse8.v v8, (t1), t2
    SET_TEST_NUM 1
    vlse8.v v16, (t1), t2
    la t1, witness_buf
    vse8.v v16, (t1)
    CHECK_MEM witness_buf, tc1_data, 4
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .byte 0x0a, 0x14, 0x1e, 0x28

.align 4
result_buf:  .space 256
witness_buf: .space 256

