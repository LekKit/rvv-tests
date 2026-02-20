/* Auto-generated test for vl4re64.v
 * Whole register load (4 reg, eew=64)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vl4re64.v: round-trip
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_data
    SAVE_CSRS
    vl4re64.v v8, (t1)
    SET_TEST_NUM 1
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 32
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 3
tc1_data:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004

.align 4
result_buf:  .space 256
witness_buf: .space 256

