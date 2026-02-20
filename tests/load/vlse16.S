/* Auto-generated test for vlse16.v
 * Strided load of 16-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlse16.v stride=4: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_data
    li t2, 4
    SAVE_CSRS
    vlse16.v v8, (t1), t2
    SET_TEST_NUM 1
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .half 0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0006, 0x0007, 0x0008
tc1_exp:
    .half 0x0001, 0x0003, 0x0005, 0x0007

.align 4
result_buf:  .space 256
witness_buf: .space 256

