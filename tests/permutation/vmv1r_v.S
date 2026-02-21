/* Auto-generated test for vmv1r.v
 * Whole register move (1 register)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmv1r.v: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    vle32.v v16, (t1)
    SAVE_CSRS
    vmv1r.v v8, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_src, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0

.align 4
result_buf:  .space 256
witness_buf: .space 256

