/* Auto-generated test for vfmv.s.f
 * Move FP scalar to vector element 0
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfmv.s.f e32: result
 *     2 = vfmv.s.f e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vmv.v.i v8, 0
    la t1, tc1_fsc
    flw fa0, 0(t1)
    SAVE_CSRS
    vfmv.s.f v8, fa0
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED_FP

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    vmv.v.i v8, 0
    la t1, tc2_fsc
    fld fa0, 0(t1)
    SAVE_CSRS
    vfmv.s.f v8, fa0
    SET_TEST_NUM 2
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_fsc:
    .word 0x402df3b6
tc1_exp:
    .word 0x402df3b6
.align 3
tc2_fsc:
    .dword 0x4005be76c8b43958
tc2_exp:
    .dword 0x4005be76c8b43958

.align 4
result_buf:  .space 256
witness_buf: .space 256

