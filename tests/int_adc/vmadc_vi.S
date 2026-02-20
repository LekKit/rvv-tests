/* Auto-generated test for vmadc.vi
 * Carry-out detection (no carry-in, vi)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmadc.vi e8 imm=0: result
 *     2 = vmadc.vi e8 imm=1: result
 *     3 = vmadc.vi e8 imm=15: result
 *     4 = vmadc.vi e8 imm=-1: result
 *     5 = vmadc.vi e16 imm=0: result
 *     6 = vmadc.vi e16 imm=1: result
 *     7 = vmadc.vi e16 imm=15: result
 *     8 = vmadc.vi e16 imm=-1: result
 *     9 = vmadc.vi e32 imm=0: result
 *    10 = vmadc.vi e32 imm=1: result
 *    11 = vmadc.vi e32 imm=15: result
 *    12 = vmadc.vi e32 imm=-1: result
 *    13 = vmadc.vi e64 imm=0: result
 *    14 = vmadc.vi e64 imm=1: result
 *    15 = vmadc.vi e64 imm=15: result
 *    16 = vmadc.vi e64 imm=-1: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 0
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
    SAVE_CSRS
    vmadc.vi v8, v16, 1
    SET_TEST_NUM 2
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 1
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 15
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
    SAVE_CSRS
    vmadc.vi v8, v16, -1
    SET_TEST_NUM 4
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 11
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 0
    SET_TEST_NUM 5
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 1
    SET_TEST_NUM 6
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 1
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 15
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
    SAVE_CSRS
    vmadc.vi v8, v16, -1
    SET_TEST_NUM 8
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 11
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 0
    SET_TEST_NUM 9
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc10_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 1
    SET_TEST_NUM 10
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 1
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 15
    SET_TEST_NUM 11
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc12_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, -1
    SET_TEST_NUM 12
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 11
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc13_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 0
    SET_TEST_NUM 13
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc14_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 1
    SET_TEST_NUM 14
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 1
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc15_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vmadc.vi v8, v16, 15
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
    SAVE_CSRS
    vmadc.vi v8, v16, -1
    SET_TEST_NUM 16
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 11
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0xff, 0xfe, 0x00, 0x01
.align 1
tc2_s2:
    .byte 0xff, 0xfe, 0x00, 0x01
.align 1
tc3_s2:
    .byte 0xff, 0xfe, 0x00, 0x01
.align 1
tc4_s2:
    .byte 0xff, 0xfe, 0x00, 0x01
.align 1
tc5_s2:
    .half 0xffff, 0xfffe, 0x0000, 0x0001
.align 1
tc6_s2:
    .half 0xffff, 0xfffe, 0x0000, 0x0001
.align 1
tc7_s2:
    .half 0xffff, 0xfffe, 0x0000, 0x0001
.align 1
tc8_s2:
    .half 0xffff, 0xfffe, 0x0000, 0x0001
.align 2
tc9_s2:
    .word 0xffffffff, 0xfffffffe, 0x00000000, 0x00000001
.align 2
tc10_s2:
    .word 0xffffffff, 0xfffffffe, 0x00000000, 0x00000001
.align 2
tc11_s2:
    .word 0xffffffff, 0xfffffffe, 0x00000000, 0x00000001
.align 2
tc12_s2:
    .word 0xffffffff, 0xfffffffe, 0x00000000, 0x00000001
.align 3
tc13_s2:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000, 0x0000000000000001
.align 3
tc14_s2:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000, 0x0000000000000001
.align 3
tc15_s2:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000, 0x0000000000000001
.align 3
tc16_s2:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000, 0x0000000000000001

.align 4
result_buf:  .space 256
witness_buf: .space 256

