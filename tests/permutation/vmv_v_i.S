/* Auto-generated test for vmv.v.i
 * Immediate broadcast to vector
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmv.v.i e8 imm=0: result
 *     2 = vmv.v.i e8 imm=1: result
 *     3 = vmv.v.i e8 imm=5: result
 *     4 = vmv.v.i e8 imm=-1: result
 *     5 = vmv.v.i e8 imm=15: result
 *     6 = vmv.v.i e8 imm=-16: result
 *     7 = vmv.v.i e16 imm=0: result
 *     8 = vmv.v.i e16 imm=1: result
 *     9 = vmv.v.i e16 imm=5: result
 *    10 = vmv.v.i e16 imm=-1: result
 *    11 = vmv.v.i e16 imm=15: result
 *    12 = vmv.v.i e16 imm=-16: result
 *    13 = vmv.v.i e32 imm=0: result
 *    14 = vmv.v.i e32 imm=1: result
 *    15 = vmv.v.i e32 imm=5: result
 *    16 = vmv.v.i e32 imm=-1: result
 *    17 = vmv.v.i e32 imm=15: result
 *    18 = vmv.v.i e32 imm=-16: result
 *    19 = vmv.v.i e64 imm=0: result
 *    20 = vmv.v.i e64 imm=1: result
 *    21 = vmv.v.i e64 imm=5: result
 *    22 = vmv.v.i e64 imm=-1: result
 *    23 = vmv.v.i e64 imm=15: result
 *    24 = vmv.v.i e64 imm=-16: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 1
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 5
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, -1
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 15
    SET_TEST_NUM 5
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, -16
    SET_TEST_NUM 6
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 0
    SET_TEST_NUM 7
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 1
    SET_TEST_NUM 8
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 5
    SET_TEST_NUM 9
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, -1
    SET_TEST_NUM 10
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 15
    SET_TEST_NUM 11
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, -16
    SET_TEST_NUM 12
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 0
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 1
    SET_TEST_NUM 14
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc14_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 5
    SET_TEST_NUM 15
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, -1
    SET_TEST_NUM 16
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 15
    SET_TEST_NUM 17
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, -16
    SET_TEST_NUM 18
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc18_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 0
    SET_TEST_NUM 19
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 1
    SET_TEST_NUM 20
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc20_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 5
    SET_TEST_NUM 21
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc21_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, -1
    SET_TEST_NUM 22
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc22_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, 15
    SET_TEST_NUM 23
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc23_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    SAVE_CSRS
    vmv.v.i v8, -16
    SET_TEST_NUM 24
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc24_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_exp:
    .byte 0x00, 0x00, 0x00, 0x00
.align 1
tc2_exp:
    .byte 0x01, 0x01, 0x01, 0x01
.align 1
tc3_exp:
    .byte 0x05, 0x05, 0x05, 0x05
.align 1
tc4_exp:
    .byte 0xff, 0xff, 0xff, 0xff
.align 1
tc5_exp:
    .byte 0x0f, 0x0f, 0x0f, 0x0f
.align 1
tc6_exp:
    .byte 0xf0, 0xf0, 0xf0, 0xf0
.align 1
tc7_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 1
tc8_exp:
    .half 0x0001, 0x0001, 0x0001, 0x0001
.align 1
tc9_exp:
    .half 0x0005, 0x0005, 0x0005, 0x0005
.align 1
tc10_exp:
    .half 0xffff, 0xffff, 0xffff, 0xffff
.align 1
tc11_exp:
    .half 0x000f, 0x000f, 0x000f, 0x000f
.align 1
tc12_exp:
    .half 0xfff0, 0xfff0, 0xfff0, 0xfff0
.align 2
tc13_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
.align 2
tc14_exp:
    .word 0x00000001, 0x00000001, 0x00000001, 0x00000001
.align 2
tc15_exp:
    .word 0x00000005, 0x00000005, 0x00000005, 0x00000005
.align 2
tc16_exp:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
.align 2
tc17_exp:
    .word 0x0000000f, 0x0000000f, 0x0000000f, 0x0000000f
.align 2
tc18_exp:
    .word 0xfffffff0, 0xfffffff0, 0xfffffff0, 0xfffffff0
.align 3
tc19_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc20_exp:
    .dword 0x0000000000000001, 0x0000000000000001, 0x0000000000000001, 0x0000000000000001
.align 3
tc21_exp:
    .dword 0x0000000000000005, 0x0000000000000005, 0x0000000000000005, 0x0000000000000005
.align 3
tc22_exp:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
.align 3
tc23_exp:
    .dword 0x000000000000000f, 0x000000000000000f, 0x000000000000000f, 0x000000000000000f
.align 3
tc24_exp:
    .dword 0xfffffffffffffff0, 0xfffffffffffffff0, 0xfffffffffffffff0, 0xfffffffffffffff0

.align 4
result_buf:  .space 256
witness_buf: .space 256

