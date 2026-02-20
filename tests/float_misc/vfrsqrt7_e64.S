/* Auto-generated test for vfrsqrt7.v
 * FP reciprocal sqrt estimate SEW=64
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfrsqrt7.v e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_src
    vle64.v v16, (t1)
    SAVE_CSRS
    vfrsqrt7.v v8, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 3
tc1_src:
    .dword 0x3ff0000000000000, 0x4010000000000000, 0x3fd0000000000000, 0x4030000000000000
tc1_exp:
    .dword 0x3fefe00000000000, 0x3fdfe00000000000, 0x3fffe00000000000, 0x3fcfe00000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

