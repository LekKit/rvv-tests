/* Auto-generated test for vsm.v
 * Unit-stride mask store
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsm.v: basic mask store
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_mask
    vlm.v v8, (t1)
    SAVE_CSRS
    la t1, result_buf
    vsm.v v8, (t1)
    SET_TEST_NUM 1
    lbu t2, 0(t1)
    andi t2, t2, 0xF
    li t3, 0b0101
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_mask:
    .byte 0b0101

.align 4
result_buf:  .space 256
witness_buf: .space 256

