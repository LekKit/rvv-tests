/* Auto-generated test for vle16ff.v
 * Fault-only-first load 16-bit
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vle16ff.v: basic (no fault)
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_data
    SAVE_CSRS
    vle16ff.v v8, (t1)
    SET_TEST_NUM 1
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 8
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .half 0x000a, 0x0014, 0x001e, 0x0028

.align 4
result_buf:  .space 256
witness_buf: .space 256

