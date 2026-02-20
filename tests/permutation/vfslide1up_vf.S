/* Auto-generated test for vfslide1up.vf
 * FP slide1 up
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfslide1up.vf e32: result
 *     2 = vfslide1up.vf e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    vle32.v v16, (t1)
    la t1, tc1_fsc
    flw fa0, 0(t1)
    SAVE_CSRS
    vfslide1up.vf v8, v16, fa0
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    CHECK_CSRS_UNCHANGED_FP

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc2_src
    vle64.v v16, (t1)
    la t1, tc2_fsc
    fld fa0, 0(t1)
    SAVE_CSRS
    vfslide1up.vf v8, v16, fa0
    SET_TEST_NUM 2
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 32
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_fsc:
    .word 0x42c60000
tc1_exp:
    .word 0x42c60000, 0x3f800000, 0x40000000, 0x40400000
.align 3
tc2_src:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc2_fsc:
    .dword 0x4058c00000000000
tc2_exp:
    .dword 0x4058c00000000000, 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

