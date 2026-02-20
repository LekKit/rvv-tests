/* Auto-generated test for vadc.vxm
 * Add with carry (vxm)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vadc.vxm e8 basic: result
 *     2 = vadc.vxm e8 carry: result
 *     3 = vadc.vxm e8 overflow: result
 *     4 = vadc.vxm e8 mixed: result
 *     5 = vadc.vxm e16 basic: result
 *     6 = vadc.vxm e16 carry: result
 *     7 = vadc.vxm e16 overflow: result
 *     8 = vadc.vxm e16 mixed: result
 *     9 = vadc.vxm e32 basic: result
 *    10 = vadc.vxm e32 carry: result
 *    11 = vadc.vxm e32 overflow: result
 *    12 = vadc.vxm e32 mixed: result
 *    13 = vadc.vxm e64 basic: result
 *    14 = vadc.vxm e64 carry: result
 *    15 = vadc.vxm e64 overflow: result
 *    16 = vadc.vxm e64 mixed: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    li a0, 0x0a
    la t1, tc1_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_s2
    vle8.v v16, (t1)
    li a0, 0x0a
    la t1, tc2_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    li a0, 0x01
    la t1, tc3_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    li a0, 0xff
    la t1, tc4_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_s2
    vle16.v v16, (t1)
    li a0, 0x000a
    la t1, tc5_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 5
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_s2
    vle16.v v16, (t1)
    li a0, 0x000a
    la t1, tc6_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 6
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_s2
    vle16.v v16, (t1)
    li a0, 0x0001
    la t1, tc7_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 7
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc8_s2
    vle16.v v16, (t1)
    li a0, 0xffff
    la t1, tc8_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 8
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_s2
    vle32.v v16, (t1)
    li a0, 0x0000000a
    la t1, tc9_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 9
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc10_s2
    vle32.v v16, (t1)
    li a0, 0x0000000a
    la t1, tc10_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 10
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_s2
    vle32.v v16, (t1)
    li a0, 0x00000001
    la t1, tc11_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 11
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc12_s2
    vle32.v v16, (t1)
    li a0, 0xffffffff
    la t1, tc12_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 12
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc13_s2
    vle64.v v16, (t1)
    li a0, 0x000000000000000a
    la t1, tc13_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 13
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc14_s2
    vle64.v v16, (t1)
    li a0, 0x000000000000000a
    la t1, tc14_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 14
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc14_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc15_s2
    vle64.v v16, (t1)
    li a0, 0x0000000000000001
    la t1, tc15_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 15
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_s2
    vle64.v v16, (t1)
    li a0, 0xffffffffffffffff
    la t1, tc16_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vxm v8, v16, a0, v0
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_mask:
    .byte 0
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_exp:
    .byte 0x0b, 0x0c, 0x0d, 0x0e
.align 1
tc2_mask:
    .byte 15
.align 1
tc2_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc2_exp:
    .byte 0x0c, 0x0d, 0x0e, 0x0f
.align 1
tc3_mask:
    .byte 10
.align 1
tc3_s2:
    .byte 0xff, 0xff, 0xfe, 0x00
tc3_exp:
    .byte 0x00, 0x01, 0xff, 0x02
.align 1
tc4_mask:
    .byte 5
.align 1
tc4_s2:
    .byte 0x00, 0xff, 0x01, 0x7f
tc4_exp:
    .byte 0x00, 0xfe, 0x01, 0x7e
.align 1
tc5_mask:
    .byte 0
.align 1
tc5_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc5_exp:
    .half 0x000b, 0x000c, 0x000d, 0x000e
.align 1
tc6_mask:
    .byte 15
.align 1
tc6_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc6_exp:
    .half 0x000c, 0x000d, 0x000e, 0x000f
.align 1
tc7_mask:
    .byte 10
.align 1
tc7_s2:
    .half 0xffff, 0xffff, 0xfffe, 0x0000
tc7_exp:
    .half 0x0000, 0x0001, 0xffff, 0x0002
.align 1
tc8_mask:
    .byte 5
.align 1
tc8_s2:
    .half 0x0000, 0xffff, 0x0001, 0x7fff
tc8_exp:
    .half 0x0000, 0xfffe, 0x0001, 0x7ffe
.align 1
tc9_mask:
    .byte 0
.align 2
tc9_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc9_exp:
    .word 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e
.align 1
tc10_mask:
    .byte 15
.align 2
tc10_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc10_exp:
    .word 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f
.align 1
tc11_mask:
    .byte 10
.align 2
tc11_s2:
    .word 0xffffffff, 0xffffffff, 0xfffffffe, 0x00000000
tc11_exp:
    .word 0x00000000, 0x00000001, 0xffffffff, 0x00000002
.align 1
tc12_mask:
    .byte 5
.align 2
tc12_s2:
    .word 0x00000000, 0xffffffff, 0x00000001, 0x7fffffff
tc12_exp:
    .word 0x00000000, 0xfffffffe, 0x00000001, 0x7ffffffe
.align 1
tc13_mask:
    .byte 0
.align 3
tc13_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc13_exp:
    .dword 0x000000000000000b, 0x000000000000000c, 0x000000000000000d, 0x000000000000000e
.align 1
tc14_mask:
    .byte 15
.align 3
tc14_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc14_exp:
    .dword 0x000000000000000c, 0x000000000000000d, 0x000000000000000e, 0x000000000000000f
.align 1
tc15_mask:
    .byte 10
.align 3
tc15_s2:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000
tc15_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0xffffffffffffffff, 0x0000000000000002
.align 1
tc16_mask:
    .byte 5
.align 3
tc16_s2:
    .dword 0x0000000000000000, 0xffffffffffffffff, 0x0000000000000001, 0x7fffffffffffffff
tc16_exp:
    .dword 0x0000000000000000, 0xfffffffffffffffe, 0x0000000000000001, 0x7ffffffffffffffe

.align 4
result_buf:  .space 256
witness_buf: .space 256

