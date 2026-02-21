/* Auto-generated test for vfnmadd.vf
 * FP FMA vfnmadd.vf
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfnmadd.vf e32 basic: wrong result
 *     2 = vfnmadd.vf e32 basic: witness changed
 *     3 = vfnmadd.vf e32 basic: CSR side-effect
 *     4 = vfnmadd.vf e32 accum: wrong result
 *     5 = vfnmadd.vf e32 accum: witness changed
 *     6 = vfnmadd.vf e32 accum: CSR side-effect
 *     7 = vfnmadd.vf e64 basic: wrong result
 *     8 = vfnmadd.vf e64 basic: witness changed
 *     9 = vfnmadd.vf e64 basic: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vfnmadd.vf SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_vd
    vle32.v v8, (t1)
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_fsc
    flw fa0, 0(t1)
    la t1, tc1_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfnmadd.vf v8, fa0, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    SET_TEST_NUM 2
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 16
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED_FP

    /* Test 4-6: vfnmadd.vf SEW=32 accum */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc4_vd
    vle32.v v8, (t1)
    la t1, tc4_s2
    vle32.v v16, (t1)
    la t1, tc4_fsc
    flw fa0, 0(t1)
    la t1, tc4_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfnmadd.vf v8, fa0, v16
    SET_TEST_NUM 4
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 16
    SET_TEST_NUM 5
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 16
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED_FP

    /* Test 7-9: vfnmadd.vf SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc7_vd
    vle64.v v8, (t1)
    la t1, tc7_s2
    vle64.v v16, (t1)
    la t1, tc7_fsc
    fld fa0, 0(t1)
    la t1, tc7_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vfnmadd.vf v8, fa0, v16
    SET_TEST_NUM 7
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 32
    SET_TEST_NUM 8
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 32
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_vd:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc1_s2:
    .word 0x41200000, 0x41a00000, 0x41f00000, 0x42200000
tc1_fsc:
    .word 0x40000000
tc1_exp:
    .word 0xc1200000, 0xc1a00000, 0xc1f00000, 0xc2200000
tc1_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc4_vd:
    .word 0x42c80000, 0x43480000, 0x43960000, 0x43c80000
tc4_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc4_fsc:
    .word 0x3f800000
tc4_exp:
    .word 0xc2ca0000, 0xc34a0000, 0xc3978000, 0xc3ca0000
tc4_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc7_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc7_s2:
    .dword 0x4024000000000000, 0x4034000000000000, 0x403e000000000000, 0x4044000000000000
tc7_fsc:
    .dword 0x4000000000000000
tc7_exp:
    .dword 0xc024000000000000, 0xc034000000000000, 0xc03e000000000000, 0xc044000000000000
tc7_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

