/* Auto-generated test for vfwmsac.vf
 * Widening FP FMA vfwmsac.vf
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfwmsac.vf basic: result
 *     2 = vfwmsac.vf basic: CSR side-effect
 *     3 = vfwmsac.vf accum: result
 *     4 = vfwmsac.vf accum: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc1_vd
    vle64.v v8, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_fsc
    flw fa0, 0(t1)
    SAVE_CSRS
    vfwmsac.vf v8, fa0, v16
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
    la t1, tc3_s2
    vle32.v v16, (t1)
    la t1, tc3_fsc
    flw fa0, 0(t1)
    SAVE_CSRS
    vfwmsac.vf v8, fa0, v16
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
tc1_s2:
    .word 0x41200000, 0x41a00000, 0x41f00000, 0x42200000
tc1_fsc:
    .word 0x40000000
tc1_exp:
    .dword 0x4034000000000000, 0x4044000000000000, 0x404e000000000000, 0x4054000000000000
.align 3
tc3_vd:
    .dword 0x4059000000000000, 0x4069000000000000, 0x4072c00000000000, 0x4079000000000000
tc3_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc3_fsc:
    .word 0x3f800000
tc3_exp:
    .dword 0xc058c00000000000, 0xc068c00000000000, 0xc072900000000000, 0xc078c00000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

