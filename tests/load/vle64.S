/* Auto-generated test for vle64.v
 * Unit-stride load of 64-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vle64.v: basic load-store round-trip
 *     2 = vle64.v: scalar verify element 0
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_data
    vle64.v v8, (t1)
    SAVE_CSRS
    SET_TEST_NUM 1
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 32
    CHECK_CSRS_UNCHANGED

    SET_TEST_NUM 2
    la t1, result_buf
    ld t2, 0(t1)
    li t3, 1
    FAIL_IF_NE t2, t3

    PASS_TEST

.data
.align 3
tc1_data:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004

.align 4
result_buf:  .space 256
witness_buf: .space 256

