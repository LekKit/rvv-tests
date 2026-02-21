/* Auto-generated test for vfsub.vf
 * FP vfsub.vf
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfsub.vf e32 basic: wrong result
 *     2 = vfsub.vf e32 basic: witness changed
 *     3 = vfsub.vf e32 basic: CSR side-effect
 *     4 = vfsub.vf e32 zero: wrong result
 *     5 = vfsub.vf e32 zero: witness changed
 *     6 = vfsub.vf e32 zero: CSR side-effect
 *     7 = vfsub.vf e32 neg: wrong result
 *     8 = vfsub.vf e32 neg: witness changed
 *     9 = vfsub.vf e32 neg: CSR side-effect
 *    10 = vfsub.vf e32 basic overlap(vd=vs2): wrong result
 *    11 = vfsub.vf e32 basic overlap(vd=vs2): witness changed
 *    12 = vfsub.vf e32 basic overlap(vd=vs2): CSR side-effect
 *    13 = vfsub.vf e32 basic masked: wrong result
 *    14 = vfsub.vf e32 basic masked: witness changed
 *    15 = vfsub.vf e32 basic masked: CSR side-effect
 *    16 = vfsub.vf e64 basic: wrong result
 *    17 = vfsub.vf e64 basic: witness changed
 *    18 = vfsub.vf e64 basic: CSR side-effect
 *    19 = vfsub.vf e64 zero: wrong result
 *    20 = vfsub.vf e64 zero: witness changed
 *    21 = vfsub.vf e64 zero: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vfsub.vf SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_fsc
    flw fa0, 0(t1)
    la t1, tc1_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfsub.vf v8, v16, fa0
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

    /* Test 4-6: vfsub.vf SEW=32 zero */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc4_s2
    vle32.v v16, (t1)
    la t1, tc4_fsc
    flw fa0, 0(t1)
    la t1, tc4_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfsub.vf v8, v16, fa0
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

    /* Test 7-9: vfsub.vf SEW=32 neg */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc7_s2
    vle32.v v16, (t1)
    la t1, tc7_fsc
    flw fa0, 0(t1)
    la t1, tc7_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfsub.vf v8, v16, fa0
    SET_TEST_NUM 7
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 16
    SET_TEST_NUM 8
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 16
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED_FP

    /* Test 10-12: vfsub.vf SEW=32 basic (vd=vs2) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc10_s2
    vle32.v v8, (t1)
    la t1, tc10_fsc
    flw fa0, 0(t1)
    la t1, tc10_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfsub.vf v8, v8, fa0
    SET_TEST_NUM 10
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 16
    SET_TEST_NUM 11
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 16
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED_FP

    /* Test 13-15: vfsub.vf SEW=32 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_s2
    vle32.v v16, (t1)
    la t1, tc13_fsc
    flw fa0, 0(t1)
    la t1, tc13_vd
    vle32.v v8, (t1)
    la t1, tc13_mask
    vlm.v v0, (t1)
    la t1, tc13_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfsub.vf v8, v16, fa0, v0.t
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16
    SET_TEST_NUM 14
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 16
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED_FP

    /* Test 16-18: vfsub.vf SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_s2
    vle64.v v16, (t1)
    la t1, tc16_fsc
    fld fa0, 0(t1)
    la t1, tc16_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vfsub.vf v8, v16, fa0
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    SET_TEST_NUM 17
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc16_w, 32
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED_FP

    /* Test 19-21: vfsub.vf SEW=64 zero */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc19_s2
    vle64.v v16, (t1)
    la t1, tc19_fsc
    fld fa0, 0(t1)
    la t1, tc19_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vfsub.vf v8, v16, fa0
    SET_TEST_NUM 19
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 32
    SET_TEST_NUM 20
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc19_w, 32
    SET_TEST_NUM 21
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_fsc:
    .word 0x41200000
tc1_exp:
    .word 0xc1100000, 0xc1000000, 0xc0e00000, 0xc0c00000
tc1_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc4_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc4_fsc:
    .word 0x00000000
tc4_exp:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc4_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc7_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc7_fsc:
    .word 0xbf800000
tc7_exp:
    .word 0x40000000, 0x40400000, 0x40800000, 0x40a00000
tc7_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc10_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc10_fsc:
    .word 0x41200000
tc10_exp:
    .word 0xc1100000, 0xc1000000, 0xc0e00000, 0xc0c00000
tc10_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 1
tc13_mask:
    .byte 10
.align 2
tc13_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc13_fsc:
    .word 0x41200000
tc13_vd:
    .word 0x42c80000, 0x42c80000, 0x42c80000, 0x42c80000
tc13_exp:
    .word 0x42c80000, 0xc1000000, 0x42c80000, 0xc0c00000
tc13_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc16_s2:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc16_fsc:
    .dword 0x4024000000000000
tc16_exp:
    .dword 0xc022000000000000, 0xc020000000000000, 0xc01c000000000000, 0xc018000000000000
tc16_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc19_s2:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc19_fsc:
    .dword 0x0000000000000000
tc19_exp:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc19_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

