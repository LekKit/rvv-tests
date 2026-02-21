/* Auto-generated test for vfwmacc.vv
 * Widening FP FMA vfwmacc.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfwmacc.vv basic: result
 *     2 = vfwmacc.vv basic: CSR side-effect
 *     3 = vfwmacc.vv accum: result
 *     4 = vfwmacc.vv accum: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc1_vd
    vle64.v v8, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s1
    vle32.v v20, (t1)
    la t1, tc1_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vfwmacc.vv v8, v20, v16
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED_FP
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc3_vd
    vle64.v v8, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_s1
    vle32.v v20, (t1)
    la t1, tc3_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vfwmacc.vv v8, v20, v16
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED_FP
    SET_TEST_NUM 3
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 32

    PASS_TEST

.data
.align 3
tc1_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc1_s1:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_s2:
    .word 0x41200000, 0x41a00000, 0x41f00000, 0x42200000
tc1_exp:
    .dword 0x4024000000000000, 0x4044000000000000, 0x4056800000000000, 0x4064000000000000
.align 3
tc3_vd:
    .dword 0x4059000000000000, 0x4069000000000000, 0x4072c00000000000, 0x4079000000000000
tc3_s1:
    .word 0x3f800000, 0x3f800000, 0x3f800000, 0x3f800000
tc3_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc3_exp:
    .dword 0x4059400000000000, 0x4069400000000000, 0x4072f00000000000, 0x4079400000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

