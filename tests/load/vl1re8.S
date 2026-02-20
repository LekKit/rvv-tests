/* Auto-generated test for vl1re8.v
 * Whole register load (1 reg, eew=8)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vl1re8.v: round-trip
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_data
    SAVE_CSRS
    vl1re8.v v8, (t1)
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 4
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 1
tc1_data:
    .byte 0x01, 0x02, 0x03, 0x04

.align 4
result_buf:  .space 256
witness_buf: .space 256

