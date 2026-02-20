/* Auto-generated test for vfmv.f.s
 * Move vector FP element 0 to FP register
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfmv.f.s e32: result
 *     2 = vfmv.f.s e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    vle32.v v16, (t1)
    SAVE_CSRS
    vfmv.f.s fa1, v16
    SET_TEST_NUM 1
    la t1, result_buf
    fsw fa1, 0(t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED_FP

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc2_src
    vle64.v v16, (t1)
    SAVE_CSRS
    vfmv.f.s fa1, v16
    SET_TEST_NUM 2
    la t1, result_buf
    fsd fa1, 0(t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0x4048f5c3, 0x00000000, 0x00000000, 0x00000000
tc1_exp:
    .word 0x4048f5c3
.align 3
tc2_src:
    .dword 0x40091eb851eb851f, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc2_exp:
    .dword 0x40091eb851eb851f

.align 4
result_buf:  .space 256
witness_buf: .space 256

