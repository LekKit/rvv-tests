/* Auto-generated test for vfmsac.vv
 * FP FMA vfmsac.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfmsac.vv e32 basic: wrong result
 *     2 = vfmsac.vv e32 basic: CSR side-effect
 *     3 = vfmsac.vv e32 accum: wrong result
 *     4 = vfmsac.vv e32 accum: CSR side-effect
 *     5 = vfmsac.vv e64 basic: wrong result
 *     6 = vfmsac.vv e64 basic: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vfmsac.vv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_vd
    vle32.v v8, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    la t1, tc1_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vfmsac.vv v8, v20, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED_FP

    /* Test 3-4: vfmsac.vv SEW=32 accum */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_vd
    vle32.v v8, (t1)
    la t1, tc3_s1
    vle32.v v20, (t1)
    la t1, tc3_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vfmsac.vv v8, v20, v16
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 16
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED_FP

    /* Test 5-6: vfmsac.vv SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc5_vd
    vle64.v v8, (t1)
    la t1, tc5_s1
    vle64.v v20, (t1)
    la t1, tc5_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vfmsac.vv v8, v20, v16
    SET_TEST_NUM 5
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 32
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_vd:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc1_s1:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_s2:
    .word 0x41200000, 0x41a00000, 0x41f00000, 0x42200000
tc1_exp:
    .word 0x41200000, 0x42200000, 0x42b40000, 0x43200000
.align 2
tc3_vd:
    .word 0x42c80000, 0x43480000, 0x43960000, 0x43c80000
tc3_s1:
    .word 0x3f800000, 0x3f800000, 0x3f800000, 0x3f800000
tc3_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc3_exp:
    .word 0xc2c60000, 0xc3460000, 0xc3948000, 0xc3c60000
.align 3
tc5_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc5_s1:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc5_s2:
    .dword 0x4024000000000000, 0x4034000000000000, 0x403e000000000000, 0x4044000000000000
tc5_exp:
    .dword 0x4024000000000000, 0x4044000000000000, 0x4056800000000000, 0x4064000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

