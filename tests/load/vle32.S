/* Auto-generated test for vle32.v
 * Unit-stride load of 32-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vle32.v: basic load-store round-trip
 *     2 = vle32.v: scalar verify element 0
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_data
    vle32.v v8, (t1)
    SAVE_CSRS
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 16
    CHECK_CSRS_UNCHANGED

    SET_TEST_NUM 2
    la t1, result_buf
    lwu t2, 0(t1)
    li t3, 1
    FAIL_IF_NE t2, t3

    PASS_TEST

.data
.align 2
tc1_data:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004

.align 4
result_buf:  .space 256
witness_buf: .space 256

