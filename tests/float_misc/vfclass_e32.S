/* Auto-generated test for vfclass.v
 * FP classify SEW=32
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfclass.v e32: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    vle32.v v16, (t1)
    SAVE_CSRS
    vfclass.v v8, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0x3f800000, 0xbf800000, 0x00000000, 0x7f800000
tc1_exp:
    .word 0x00000040, 0x00000002, 0x00000010, 0x00000080

.align 4
result_buf:  .space 256
witness_buf: .space 256

