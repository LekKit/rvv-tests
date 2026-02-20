/* Auto-generated test for vfwcvt.rtz.x.f.v
 * Widening FP conversion vfwcvt.rtz.x.f.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfwcvt.rtz.x.f.v: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    vle32.v v16, (t1)
    SAVE_CSRS
    vfwcvt.rtz.x.f.v v8, v16
    SET_TEST_NUM 1
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32
    CHECK_VSTART_ZERO
    CHECK_FCSR_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_src:
    .word 0x00000000, 0x3f800000, 0x42280000, 0x42c80000
tc1_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0x000000000000002a, 0x0000000000000064

.align 4
result_buf:  .space 256
witness_buf: .space 256

