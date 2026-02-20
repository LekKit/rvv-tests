/* Auto-generated test for vfmv.v.f
 * FP scalar broadcast to vector
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfmv.v.f e32: result
 *     2 = vfmv.v.f e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_fsc
    flw fa0, 0(t1)
    SAVE_CSRS
    vfmv.v.f v8, fa0
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    CHECK_CSRS_UNCHANGED_FP

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc2_fsc
    fld fa0, 0(t1)
    SAVE_CSRS
    vfmv.v.f v8, fa0
    SET_TEST_NUM 2
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 32
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_fsc:
    .word 0x40f00000
tc1_exp:
    .word 0x40f00000, 0x40f00000, 0x40f00000, 0x40f00000
.align 3
tc2_fsc:
    .dword 0x401e000000000000
tc2_exp:
    .dword 0x401e000000000000, 0x401e000000000000, 0x401e000000000000, 0x401e000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

