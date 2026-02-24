/* Auto-generated test for vmfeq.vf
 * FP compare vmfeq.vf
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmfeq.vf e32 basic: wrong mask result
 *     2 = vmfeq.vf e32 basic: witness changed
 *     3 = vmfeq.vf e32 basic: CSR side-effect
 *     4 = vmfeq.vf e32 zero: wrong mask result
 *     5 = vmfeq.vf e32 zero: witness changed
 *     6 = vmfeq.vf e32 zero: CSR side-effect
 *     7 = vmfeq.vf e32 neg: wrong mask result
 *     8 = vmfeq.vf e32 neg: witness changed
 *     9 = vmfeq.vf e32 neg: CSR side-effect
 *    10 = vmfeq.vf e32 mixed masked: wrong mask result
 *    11 = vmfeq.vf e32 mixed masked: witness changed
 *    12 = vmfeq.vf e32 mixed masked: CSR side-effect
 *    13 = vmfeq.vf e64 basic: wrong mask result
 *    14 = vmfeq.vf e64 basic: witness changed
 *    15 = vmfeq.vf e64 basic: CSR side-effect
 *    16 = vmfeq.vf e64 zero: wrong mask result
 *    17 = vmfeq.vf e64 zero: witness changed
 *    18 = vmfeq.vf e64 zero: CSR side-effect
 *    19 = vmfeq.vf e64 mixed masked: wrong mask result
 *    20 = vmfeq.vf e64 mixed masked: witness changed
 *    21 = vmfeq.vf e64 mixed masked: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vmfeq.vf SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_fsc
    flw fa0, 0(t1)
    la t1, tc1_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmfeq.vf v8, v16, fa0
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 2
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 16
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED_FP

    /* Test 4-6: vmfeq.vf SEW=32 zero */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc4_s2
    vle32.v v16, (t1)
    la t1, tc4_fsc
    flw fa0, 0(t1)
    la t1, tc4_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmfeq.vf v8, v16, fa0
    SET_TEST_NUM 4
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 5
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 16
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED_FP

    /* Test 7-9: vmfeq.vf SEW=32 neg */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc7_s2
    vle32.v v16, (t1)
    la t1, tc7_fsc
    flw fa0, 0(t1)
    la t1, tc7_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmfeq.vf v8, v16, fa0
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 8
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 16
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED_FP

    /* Test 10-12: vmfeq.vf SEW=32 mixed (masked) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc10_s2
    vle32.v v16, (t1)
    la t1, tc10_sc
    flw fa0, 0(t1)
    la t1, tc10_vdinit
    vlm.v v8, (t1)
    la t1, tc10_mask
    vlm.v v0, (t1)
    la t1, tc10_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmfeq.vf v8, v16, fa0, v0.t
    SET_TEST_NUM 10
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 10
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 11
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 16
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED_FP

    /* Test 13-15: vmfeq.vf SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc13_s2
    vle64.v v16, (t1)
    la t1, tc13_fsc
    fld fa0, 0(t1)
    la t1, tc13_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vmfeq.vf v8, v16, fa0
    SET_TEST_NUM 13
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 14
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 32
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED_FP

    /* Test 16-18: vmfeq.vf SEW=64 zero */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_s2
    vle64.v v16, (t1)
    la t1, tc16_fsc
    fld fa0, 0(t1)
    la t1, tc16_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vmfeq.vf v8, v16, fa0
    SET_TEST_NUM 16
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 17
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc16_w, 32
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED_FP

    /* Test 19-21: vmfeq.vf SEW=64 mixed (masked) */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc19_s2
    vle64.v v16, (t1)
    la t1, tc19_sc
    fld fa0, 0(t1)
    la t1, tc19_vdinit
    vlm.v v8, (t1)
    la t1, tc19_mask
    vlm.v v0, (t1)
    la t1, tc19_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vmfeq.vf v8, v16, fa0, v0.t
    SET_TEST_NUM 19
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 10
    FAIL_IF_NE t2, t3
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
tc1_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc4_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc4_fsc:
    .word 0x00000000
tc4_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc7_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc7_fsc:
    .word 0xbf800000
tc7_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 1
tc10_mask:
    .byte 5
tc10_vdinit:
    .byte 10
.align 2
tc10_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc10_sc:
    .word 0x40200000
tc10_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc13_s2:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc13_fsc:
    .dword 0x4024000000000000
tc13_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc16_s2:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc16_fsc:
    .dword 0x0000000000000000
tc16_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 1
tc19_mask:
    .byte 5
tc19_vdinit:
    .byte 10
.align 3
tc19_s2:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc19_sc:
    .dword 0x4004000000000000
tc19_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

