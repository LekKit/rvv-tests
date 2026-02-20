/* Auto-generated test for vadc.vim
 * Add with carry (vim)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vadc.vim e8 imm=0: result
 *     2 = vadc.vim e8 imm=1: result
 *     3 = vadc.vim e8 imm=5: result
 *     4 = vadc.vim e8 imm=-1: result
 *     5 = vadc.vim e8 imm=-16: result
 *     6 = vadc.vim e16 imm=0: result
 *     7 = vadc.vim e16 imm=1: result
 *     8 = vadc.vim e16 imm=5: result
 *     9 = vadc.vim e16 imm=-1: result
 *    10 = vadc.vim e16 imm=-16: result
 *    11 = vadc.vim e32 imm=0: result
 *    12 = vadc.vim e32 imm=1: result
 *    13 = vadc.vim e32 imm=5: result
 *    14 = vadc.vim e32 imm=-1: result
 *    15 = vadc.vim e32 imm=-16: result
 *    16 = vadc.vim e64 imm=0: result
 *    17 = vadc.vim e64 imm=1: result
 *    18 = vadc.vim e64 imm=5: result
 *    19 = vadc.vim e64 imm=-1: result
 *    20 = vadc.vim e64 imm=-16: result
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
    vadc.vim v8, v16, 0, v0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_s2
    vle8.v v16, (t1)
    la t1, tc2_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 1, v0
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    la t1, tc3_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 5, v0
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    la t1, tc4_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, -1, v0
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_s2
    vle8.v v16, (t1)
    la t1, tc5_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, -16, v0
    SET_TEST_NUM 5
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_s2
    vle16.v v16, (t1)
    la t1, tc6_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 0, v0
    SET_TEST_NUM 6
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_s2
    vle16.v v16, (t1)
    la t1, tc7_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 1, v0
    SET_TEST_NUM 7
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc8_s2
    vle16.v v16, (t1)
    la t1, tc8_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 5, v0
    SET_TEST_NUM 8
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_s2
    vle16.v v16, (t1)
    la t1, tc9_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, -1, v0
    SET_TEST_NUM 9
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc10_s2
    vle16.v v16, (t1)
    la t1, tc10_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, -16, v0
    SET_TEST_NUM 10
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_s2
    vle32.v v16, (t1)
    la t1, tc11_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 0, v0
    SET_TEST_NUM 11
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc12_s2
    vle32.v v16, (t1)
    la t1, tc12_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 1, v0
    SET_TEST_NUM 12
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_s2
    vle32.v v16, (t1)
    la t1, tc13_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 5, v0
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc14_s2
    vle32.v v16, (t1)
    la t1, tc14_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, -1, v0
    SET_TEST_NUM 14
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc14_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc15_s2
    vle32.v v16, (t1)
    la t1, tc15_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, -16, v0
    SET_TEST_NUM 15
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_s2
    vle64.v v16, (t1)
    la t1, tc16_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 0, v0
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc17_s2
    vle64.v v16, (t1)
    la t1, tc17_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 1, v0
    SET_TEST_NUM 17
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc18_s2
    vle64.v v16, (t1)
    la t1, tc18_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, 5, v0
    SET_TEST_NUM 18
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc18_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc19_s2
    vle64.v v16, (t1)
    la t1, tc19_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, -1, v0
    SET_TEST_NUM 19
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc20_s2
    vle64.v v16, (t1)
    la t1, tc20_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vadc.vim v8, v16, -16, v0
    SET_TEST_NUM 20
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc20_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_mask:
    .byte 10
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0xff, 0x00
tc1_exp:
    .byte 0x01, 0x03, 0xff, 0x01
.align 1
tc2_mask:
    .byte 10
.align 1
tc2_s2:
    .byte 0x01, 0x02, 0xff, 0x00
tc2_exp:
    .byte 0x02, 0x04, 0x00, 0x02
.align 1
tc3_mask:
    .byte 10
.align 1
tc3_s2:
    .byte 0x01, 0x02, 0xff, 0x00
tc3_exp:
    .byte 0x06, 0x08, 0x04, 0x06
.align 1
tc4_mask:
    .byte 10
.align 1
tc4_s2:
    .byte 0x01, 0x02, 0xff, 0x00
tc4_exp:
    .byte 0x00, 0x02, 0xfe, 0x00
.align 1
tc5_mask:
    .byte 10
.align 1
tc5_s2:
    .byte 0x01, 0x02, 0xff, 0x00
tc5_exp:
    .byte 0xf1, 0xf3, 0xef, 0xf1
.align 1
tc6_mask:
    .byte 10
.align 1
tc6_s2:
    .half 0x0001, 0x0002, 0xffff, 0x0000
tc6_exp:
    .half 0x0001, 0x0003, 0xffff, 0x0001
.align 1
tc7_mask:
    .byte 10
.align 1
tc7_s2:
    .half 0x0001, 0x0002, 0xffff, 0x0000
tc7_exp:
    .half 0x0002, 0x0004, 0x0000, 0x0002
.align 1
tc8_mask:
    .byte 10
.align 1
tc8_s2:
    .half 0x0001, 0x0002, 0xffff, 0x0000
tc8_exp:
    .half 0x0006, 0x0008, 0x0004, 0x0006
.align 1
tc9_mask:
    .byte 10
.align 1
tc9_s2:
    .half 0x0001, 0x0002, 0xffff, 0x0000
tc9_exp:
    .half 0x0000, 0x0002, 0xfffe, 0x0000
.align 1
tc10_mask:
    .byte 10
.align 1
tc10_s2:
    .half 0x0001, 0x0002, 0xffff, 0x0000
tc10_exp:
    .half 0xfff1, 0xfff3, 0xffef, 0xfff1
.align 1
tc11_mask:
    .byte 10
.align 2
tc11_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x00000000
tc11_exp:
    .word 0x00000001, 0x00000003, 0xffffffff, 0x00000001
.align 1
tc12_mask:
    .byte 10
.align 2
tc12_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x00000000
tc12_exp:
    .word 0x00000002, 0x00000004, 0x00000000, 0x00000002
.align 1
tc13_mask:
    .byte 10
.align 2
tc13_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x00000000
tc13_exp:
    .word 0x00000006, 0x00000008, 0x00000004, 0x00000006
.align 1
tc14_mask:
    .byte 10
.align 2
tc14_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x00000000
tc14_exp:
    .word 0x00000000, 0x00000002, 0xfffffffe, 0x00000000
.align 1
tc15_mask:
    .byte 10
.align 2
tc15_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x00000000
tc15_exp:
    .word 0xfffffff1, 0xfffffff3, 0xffffffef, 0xfffffff1
.align 1
tc16_mask:
    .byte 10
.align 3
tc16_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0xffffffffffffffff, 0x0000000000000000
tc16_exp:
    .dword 0x0000000000000001, 0x0000000000000003, 0xffffffffffffffff, 0x0000000000000001
.align 1
tc17_mask:
    .byte 10
.align 3
tc17_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0xffffffffffffffff, 0x0000000000000000
tc17_exp:
    .dword 0x0000000000000002, 0x0000000000000004, 0x0000000000000000, 0x0000000000000002
.align 1
tc18_mask:
    .byte 10
.align 3
tc18_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0xffffffffffffffff, 0x0000000000000000
tc18_exp:
    .dword 0x0000000000000006, 0x0000000000000008, 0x0000000000000004, 0x0000000000000006
.align 1
tc19_mask:
    .byte 10
.align 3
tc19_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0xffffffffffffffff, 0x0000000000000000
tc19_exp:
    .dword 0x0000000000000000, 0x0000000000000002, 0xfffffffffffffffe, 0x0000000000000000
.align 1
tc20_mask:
    .byte 10
.align 3
tc20_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0xffffffffffffffff, 0x0000000000000000
tc20_exp:
    .dword 0xfffffffffffffff1, 0xfffffffffffffff3, 0xffffffffffffffef, 0xfffffffffffffff1

.align 4
result_buf:  .space 256
witness_buf: .space 256

