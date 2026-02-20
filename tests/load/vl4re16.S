/* Auto-generated test for vl4re16.v
 * Whole register load (4 reg, eew=16)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vl4re16.v: round-trip
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_data
    SAVE_CSRS
    vl4re16.v v8, (t1)
    SET_TEST_NUM 1
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 8
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 1
tc1_data:
    .half 0x0001, 0x0002, 0x0003, 0x0004

.align 4
result_buf:  .space 256
witness_buf: .space 256

