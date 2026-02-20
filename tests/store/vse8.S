/* Auto-generated test for vse8.v
 * Unit-stride store of 8-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vse8.v: store and scalar-verify
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_data
    vle8.v v8, (t1)
    SAVE_CSRS
    la t1, result_buf
    vse8.v v8, (t1)
    SET_TEST_NUM 1
    lbu t2, 0(t1)
    li t3, 1
    FAIL_IF_NE t2, t3
    lbu t2, 1(t1)
    li t3, 2
    FAIL_IF_NE t2, t3
    lbu t2, 2(t1)
    li t3, 3
    FAIL_IF_NE t2, t3
    lbu t2, 3(t1)
    li t3, 4
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .byte 0x01, 0x02, 0x03, 0x04

.align 4
result_buf:  .space 256
witness_buf: .space 256

