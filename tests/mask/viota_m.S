/* Auto-generated test for viota.m
 * Vector iota (prefix sum of mask)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = viota.m e8 all: result
 *     2 = viota.m e8 none: result
 *     3 = viota.m e8 alt: result
 *     4 = viota.m e8 first: result
 *     5 = viota.m e16 all: result
 *     6 = viota.m e16 none: result
 *     7 = viota.m e16 alt: result
 *     8 = viota.m e16 first: result
 *     9 = viota.m e32 all: result
 *    10 = viota.m e32 none: result
 *    11 = viota.m e32 alt: result
 *    12 = viota.m e32 first: result
 *    13 = viota.m e64 all: result
 *    14 = viota.m e64 none: result
 *    15 = viota.m e64 alt: result
 *    16 = viota.m e64 first: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 5
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 6
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 7
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc8_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 8
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 9
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc10_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 10
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 11
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc12_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 12
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc13_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 13
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc14_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 14
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc14_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc15_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 15
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    viota.m v8, v0
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_mask:
    .byte 15
.align 1
tc1_exp:
    .byte 0x00, 0x01, 0x02, 0x03
.align 1
tc2_mask:
    .byte 0
.align 1
tc2_exp:
    .byte 0x00, 0x00, 0x00, 0x00
.align 1
tc3_mask:
    .byte 10
.align 1
tc3_exp:
    .byte 0x00, 0x00, 0x01, 0x01
.align 1
tc4_mask:
    .byte 1
.align 1
tc4_exp:
    .byte 0x00, 0x01, 0x01, 0x01
.align 1
tc5_mask:
    .byte 15
.align 1
tc5_exp:
    .half 0x0000, 0x0001, 0x0002, 0x0003
.align 1
tc6_mask:
    .byte 0
.align 1
tc6_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 1
tc7_mask:
    .byte 10
.align 1
tc7_exp:
    .half 0x0000, 0x0000, 0x0001, 0x0001
.align 1
tc8_mask:
    .byte 1
.align 1
tc8_exp:
    .half 0x0000, 0x0001, 0x0001, 0x0001
.align 1
tc9_mask:
    .byte 15
.align 2
tc9_exp:
    .word 0x00000000, 0x00000001, 0x00000002, 0x00000003
.align 1
tc10_mask:
    .byte 0
.align 2
tc10_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
.align 1
tc11_mask:
    .byte 10
.align 2
tc11_exp:
    .word 0x00000000, 0x00000000, 0x00000001, 0x00000001
.align 1
tc12_mask:
    .byte 1
.align 2
tc12_exp:
    .word 0x00000000, 0x00000001, 0x00000001, 0x00000001
.align 1
tc13_mask:
    .byte 15
.align 3
tc13_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0x0000000000000002, 0x0000000000000003
.align 1
tc14_mask:
    .byte 0
.align 3
tc14_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 1
tc15_mask:
    .byte 10
.align 3
tc15_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000001, 0x0000000000000001
.align 1
tc16_mask:
    .byte 1
.align 3
tc16_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0x0000000000000001, 0x0000000000000001

.align 4
result_buf:  .space 256
witness_buf: .space 256

