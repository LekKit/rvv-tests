/* Auto-generated test for vfncvt.f.f.w
 * Narrowing FP conversion vfncvt.f.f.w
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfncvt.f.f.w: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc1_src
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vfncvt.f.f.w v8, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    CHECK_VSTART_ZERO
    CHECK_FCSR_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_src:
    .dword 0x0000000000000000, 0x3ff0000000000000, 0x4045000000000000, 0x4059000000000000
tc1_exp:
    .word 0x00000000, 0x3f800000, 0x42280000, 0x42c80000

.align 4
result_buf:  .space 256
witness_buf: .space 256

