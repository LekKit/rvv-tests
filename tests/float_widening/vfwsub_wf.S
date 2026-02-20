/* Auto-generated test for vfwsub.wf
 * Widening FP vfwsub.wf
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfwsub.wf basic: result
 *     2 = vfwsub.wf neg: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc1_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_fsc
    flw fa0, 0(t1)
    SAVE_CSRS
    vfwsub.wf v8, v16, fa0
    SET_TEST_NUM 1
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32
    CHECK_VSTART_ZERO
    CHECK_FCSR_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc2_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc2_fsc
    flw fa0, 0(t1)
    SAVE_CSRS
    vfwsub.wf v8, v16, fa0
    SET_TEST_NUM 2
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 32
    CHECK_VSTART_ZERO
    CHECK_FCSR_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_s2w:
    .dword 0x4024000000000000, 0x4034000000000000, 0x403e000000000000, 0x4044000000000000
tc1_fsc:
    .word 0x3f800000
tc1_exp:
    .dword 0x4022000000000000, 0x4033000000000000, 0x403d000000000000, 0x4043800000000000
.align 3
tc2_s2w:
    .dword 0xbff0000000000000, 0xc000000000000000, 0x4008000000000000, 0x4010000000000000
tc2_fsc:
    .word 0xc0400000
tc2_exp:
    .dword 0x4000000000000000, 0x3ff0000000000000, 0x4018000000000000, 0x401c000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

