/* Auto-generated test for vle32ff.v
 * Fault-only-first load 32-bit
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vle32ff.v: basic (no fault)
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_data
    SAVE_CSRS
    vle32ff.v v8, (t1)
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 16
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 2
tc1_data:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028

.align 4
result_buf:  .space 256
witness_buf: .space 256

