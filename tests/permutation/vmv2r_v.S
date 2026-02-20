/* Auto-generated test for vmv2r.v
 * Whole register move (2 registers)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmv2r.v: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 8
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc1_src
    vle32.v v16, (t1)
    SAVE_CSRS
    vmv2r.v v8, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_src, 32
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0xa0000000, 0xb1111111, 0xc2222222, 0xd3333333, 0xe4444444, 0xf5555555, 0x06666666, 0x17777777

.align 4
result_buf:  .space 256
witness_buf: .space 256

