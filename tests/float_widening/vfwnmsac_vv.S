/* Auto-generated test for vfwnmsac.vv
 * Widening FP FMA vfwnmsac.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfwnmsac.vv basic: result
 *     2 = vfwnmsac.vv accum: result
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
    vfwnmsac.vv v8, v20, v16
    SET_TEST_NUM 1
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32
    CHECK_VSTART_ZERO
    CHECK_FCSR_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc2_vd
    vle64.v v8, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc2_s1
    vle32.v v20, (t1)
    la t1, tc2_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vfwnmsac.vv v8, v20, v16
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
tc1_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc1_s1:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_s2:
    .word 0x41200000, 0x41a00000, 0x41f00000, 0x42200000
tc1_exp:
    .dword 0xc024000000000000, 0xc044000000000000, 0xc056800000000000, 0xc064000000000000
.align 3
tc2_vd:
    .dword 0x4059000000000000, 0x4069000000000000, 0x4072c00000000000, 0x4079000000000000
tc2_s1:
    .word 0x3f800000, 0x3f800000, 0x3f800000, 0x3f800000
tc2_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc2_exp:
    .dword 0x4058c00000000000, 0x4068c00000000000, 0x4072900000000000, 0x4078c00000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

