/* Auto-generated test for vl1re32.v
 * Whole register load (1 reg, eew=32)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vl1re32.v: round-trip
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_data
    SAVE_CSRS
    vl1re32.v v8, (t1)
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_data:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004

.align 4
result_buf:  .space 256
witness_buf: .space 256

