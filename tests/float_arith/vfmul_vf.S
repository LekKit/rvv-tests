/* Auto-generated test for vfmul.vf
 * FP vfmul.vf
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfmul.vf e32 basic: wrong result
 *     2 = vfmul.vf e32 basic: CSR side-effect
 *     3 = vfmul.vf e32 zero: wrong result
 *     4 = vfmul.vf e32 zero: CSR side-effect
 *     5 = vfmul.vf e32 neg: wrong result
 *     6 = vfmul.vf e32 neg: CSR side-effect
 *     7 = vfmul.vf e64 basic: wrong result
 *     8 = vfmul.vf e64 basic: CSR side-effect
 *     9 = vfmul.vf e64 zero: wrong result
 *    10 = vfmul.vf e64 zero: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vfmul.vf SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_fsc
    flw fa0, 0(t1)
    SAVE_CSRS
    vfmul.vf v8, v16, fa0
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED_FP

    /* Test 3-4: vfmul.vf SEW=32 zero */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_s2
    vle32.v v16, (t1)
    la t1, tc3_fsc
    flw fa0, 0(t1)
    SAVE_CSRS
    vfmul.vf v8, v16, fa0
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 16
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED_FP

    /* Test 5-6: vfmul.vf SEW=32 neg */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc5_s2
    vle32.v v16, (t1)
    la t1, tc5_fsc
    flw fa0, 0(t1)
    SAVE_CSRS
    vfmul.vf v8, v16, fa0
    SET_TEST_NUM 5
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 16
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED_FP

    /* Test 7-8: vfmul.vf SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc7_s2
    vle64.v v16, (t1)
    la t1, tc7_fsc
    fld fa0, 0(t1)
    SAVE_CSRS
    vfmul.vf v8, v16, fa0
    SET_TEST_NUM 7
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 32
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED_FP

    /* Test 9-10: vfmul.vf SEW=64 zero */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc9_s2
    vle64.v v16, (t1)
    la t1, tc9_fsc
    fld fa0, 0(t1)
    SAVE_CSRS
    vfmul.vf v8, v16, fa0
    SET_TEST_NUM 9
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 32
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_fsc:
    .word 0x41200000
tc1_exp:
    .word 0x41200000, 0x41a00000, 0x41f00000, 0x42200000
.align 2
tc3_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc3_fsc:
    .word 0x00000000
tc3_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
.align 2
tc5_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc5_fsc:
    .word 0xbf800000
tc5_exp:
    .word 0xbf800000, 0xc0000000, 0xc0400000, 0xc0800000
.align 3
tc7_s2:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc7_fsc:
    .dword 0x4024000000000000
tc7_exp:
    .dword 0x4024000000000000, 0x4034000000000000, 0x403e000000000000, 0x4044000000000000
.align 3
tc9_s2:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc9_fsc:
    .dword 0x0000000000000000
tc9_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

