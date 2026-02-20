/* Auto-generated test for vfwmul.vv
 * Widening FP vfwmul.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfwmul.vv basic: result
 *     2 = vfwmul.vv negpos: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vfwmul.vv v8, v16, v20
    SET_TEST_NUM 1
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32
    CHECK_VSTART_ZERO
    CHECK_FCSR_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc2_s2
    vle32.v v16, (t1)
    la t1, tc2_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vfwmul.vv v8, v16, v20
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
tc1_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_s1:
    .word 0x3f000000, 0x3e800000, 0x3e000000, 0x3d800000
tc1_exp:
    .dword 0x3fe0000000000000, 0x3fe0000000000000, 0x3fd8000000000000, 0x3fd0000000000000
.align 3
tc2_s2:
    .word 0xbf800000, 0xc0000000, 0x40400000, 0x40800000
tc2_s1:
    .word 0x3f800000, 0x40000000, 0xc0400000, 0xc0800000
tc2_exp:
    .dword 0xbff0000000000000, 0xc010000000000000, 0xc022000000000000, 0xc030000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

