/* Auto-generated test for vs4r.v
 * Whole register store (4 reg)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vs4r.v: round-trip
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_data
    vle32.v v8, (t1)
    SAVE_CSRS
    la t1, result_buf
    vs4r.v v8, (t1)
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_data, 16
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 2
tc1_data:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0

.align 4
result_buf:  .space 256
witness_buf: .space 256

