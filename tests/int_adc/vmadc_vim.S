/* Auto-generated test for vmadc.vim
 * Carry-out detection with carry-in (vim)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmadc.vim e8 imm=0: result
 *     2 = vmadc.vim e8 imm=1: result
 *     3 = vmadc.vim e8 imm=15: result
 *     4 = vmadc.vim e8 imm=-1: result
 *     5 = vmadc.vim e8 imm=-16: result
 *     6 = vmadc.vim e16 imm=0: result
 *     7 = vmadc.vim e16 imm=1: result
 *     8 = vmadc.vim e16 imm=15: result
 *     9 = vmadc.vim e16 imm=-1: result
 *    10 = vmadc.vim e16 imm=-16: result
 *    11 = vmadc.vim e32 imm=0: result
 *    12 = vmadc.vim e32 imm=1: result
 *    13 = vmadc.vim e32 imm=15: result
 *    14 = vmadc.vim e32 imm=-1: result
 *    15 = vmadc.vim e32 imm=-16: result
 *    16 = vmadc.vim e64 imm=0: result
 *    17 = vmadc.vim e64 imm=1: result
 *    18 = vmadc.vim e64 imm=15: result
 *    19 = vmadc.vim e64 imm=-1: result
 *    20 = vmadc.vim e64 imm=-16: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 0, v0
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_s2
    vle8.v v16, (t1)
    la t1, tc2_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 1, v0
    SET_TEST_NUM 2
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    la t1, tc3_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 15, v0
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    la t1, tc4_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, -1, v0
    SET_TEST_NUM 4
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 11
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_s2
    vle8.v v16, (t1)
    la t1, tc5_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, -16, v0
    SET_TEST_NUM 5
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_s2
    vle16.v v16, (t1)
    la t1, tc6_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 0, v0
    SET_TEST_NUM 6
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_s2
    vle16.v v16, (t1)
    la t1, tc7_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 1, v0
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc8_s2
    vle16.v v16, (t1)
    la t1, tc8_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 15, v0
    SET_TEST_NUM 8
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_s2
    vle16.v v16, (t1)
    la t1, tc9_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, -1, v0
    SET_TEST_NUM 9
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 11
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc10_s2
    vle16.v v16, (t1)
    la t1, tc10_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, -16, v0
    SET_TEST_NUM 10
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_s2
    vle32.v v16, (t1)
    la t1, tc11_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 0, v0
    SET_TEST_NUM 11
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc12_s2
    vle32.v v16, (t1)
    la t1, tc12_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 1, v0
    SET_TEST_NUM 12
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_s2
    vle32.v v16, (t1)
    la t1, tc13_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 15, v0
    SET_TEST_NUM 13
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc14_s2
    vle32.v v16, (t1)
    la t1, tc14_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, -1, v0
    SET_TEST_NUM 14
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 11
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc15_s2
    vle32.v v16, (t1)
    la t1, tc15_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, -16, v0
    SET_TEST_NUM 15
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_s2
    vle64.v v16, (t1)
    la t1, tc16_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 0, v0
    SET_TEST_NUM 16
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc17_s2
    vle64.v v16, (t1)
    la t1, tc17_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 1, v0
    SET_TEST_NUM 17
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc18_s2
    vle64.v v16, (t1)
    la t1, tc18_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, 15, v0
    SET_TEST_NUM 18
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc19_s2
    vle64.v v16, (t1)
    la t1, tc19_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, -1, v0
    SET_TEST_NUM 19
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 11
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc20_s2
    vle64.v v16, (t1)
    la t1, tc20_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmadc.vim v8, v16, -16, v0
    SET_TEST_NUM 20
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_mask:
    .byte 10
.align 1
tc1_s2:
    .byte 0xff, 0xfe, 0x00, 0x01
.align 1
tc2_mask:
    .byte 10
.align 1
tc2_s2:
    .byte 0xff, 0xfe, 0x00, 0x01
.align 1
tc3_mask:
    .byte 10
.align 1
tc3_s2:
    .byte 0xff, 0xfe, 0x00, 0x01
.align 1
tc4_mask:
    .byte 10
.align 1
tc4_s2:
    .byte 0xff, 0xfe, 0x00, 0x01
.align 1
tc5_mask:
    .byte 10
.align 1
tc5_s2:
    .byte 0xff, 0xfe, 0x00, 0x01
.align 1
tc6_mask:
    .byte 10
.align 1
tc6_s2:
    .half 0xffff, 0xfffe, 0x0000, 0x0001
.align 1
tc7_mask:
    .byte 10
.align 1
tc7_s2:
    .half 0xffff, 0xfffe, 0x0000, 0x0001
.align 1
tc8_mask:
    .byte 10
.align 1
tc8_s2:
    .half 0xffff, 0xfffe, 0x0000, 0x0001
.align 1
tc9_mask:
    .byte 10
.align 1
tc9_s2:
    .half 0xffff, 0xfffe, 0x0000, 0x0001
.align 1
tc10_mask:
    .byte 10
.align 1
tc10_s2:
    .half 0xffff, 0xfffe, 0x0000, 0x0001
.align 1
tc11_mask:
    .byte 10
.align 2
tc11_s2:
    .word 0xffffffff, 0xfffffffe, 0x00000000, 0x00000001
.align 1
tc12_mask:
    .byte 10
.align 2
tc12_s2:
    .word 0xffffffff, 0xfffffffe, 0x00000000, 0x00000001
.align 1
tc13_mask:
    .byte 10
.align 2
tc13_s2:
    .word 0xffffffff, 0xfffffffe, 0x00000000, 0x00000001
.align 1
tc14_mask:
    .byte 10
.align 2
tc14_s2:
    .word 0xffffffff, 0xfffffffe, 0x00000000, 0x00000001
.align 1
tc15_mask:
    .byte 10
.align 2
tc15_s2:
    .word 0xffffffff, 0xfffffffe, 0x00000000, 0x00000001
.align 1
tc16_mask:
    .byte 10
.align 3
tc16_s2:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000, 0x0000000000000001
.align 1
tc17_mask:
    .byte 10
.align 3
tc17_s2:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000, 0x0000000000000001
.align 1
tc18_mask:
    .byte 10
.align 3
tc18_s2:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000, 0x0000000000000001
.align 1
tc19_mask:
    .byte 10
.align 3
tc19_s2:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000, 0x0000000000000001
.align 1
tc20_mask:
    .byte 10
.align 3
tc20_s2:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000, 0x0000000000000001

.align 4
result_buf:  .space 256
witness_buf: .space 256

