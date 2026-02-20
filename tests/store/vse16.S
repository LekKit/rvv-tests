/* Auto-generated test for vse16.v
 * Unit-stride store of 16-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vse16.v: store and scalar-verify
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_data
    vle16.v v8, (t1)
    SAVE_CSRS
    la t1, result_buf
    vse16.v v8, (t1)
    SET_TEST_NUM 1
    lhu t2, 0(t1)
    li t3, 1
    FAIL_IF_NE t2, t3
    lhu t2, 2(t1)
    li t3, 2
    FAIL_IF_NE t2, t3
    lhu t2, 4(t1)
    li t3, 3
    FAIL_IF_NE t2, t3
    lhu t2, 6(t1)
    li t3, 4
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .half 0x0001, 0x0002, 0x0003, 0x0004

.align 4
result_buf:  .space 256
witness_buf: .space 256

