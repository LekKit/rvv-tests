/* Auto-generated test for vsadd.vi
 * Saturating vsadd.vi
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsadd.vi e8 imm=0: result
 *     2 = vsadd.vi e8 imm=1: result
 *     3 = vsadd.vi e8 imm=15: result
 *     4 = vsadd.vi e8 imm=-1: result
 *     5 = vsadd.vi e8 imm=-16: result
 *     6 = vsadd.vi e16 imm=0: result
 *     7 = vsadd.vi e16 imm=1: result
 *     8 = vsadd.vi e16 imm=15: result
 *     9 = vsadd.vi e16 imm=-1: result
 *    10 = vsadd.vi e16 imm=-16: result
 *    11 = vsadd.vi e32 imm=0: result
 *    12 = vsadd.vi e32 imm=1: result
 *    13 = vsadd.vi e32 imm=15: result
 *    14 = vsadd.vi e32 imm=-1: result
 *    15 = vsadd.vi e32 imm=-16: result
 *    16 = vsadd.vi e64 imm=0: result
 *    17 = vsadd.vi e64 imm=1: result
 *    18 = vsadd.vi e64 imm=15: result
 *    19 = vsadd.vi e64 imm=-1: result
 *    20 = vsadd.vi e64 imm=-16: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 1
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 15
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, -1
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, -16
    SET_TEST_NUM 5
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 4
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 0
    SET_TEST_NUM 6
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 8
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 1
    SET_TEST_NUM 7
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc8_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 15
    SET_TEST_NUM 8
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 8
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, -1
    SET_TEST_NUM 9
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 8
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc10_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, -16
    SET_TEST_NUM 10
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 8
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 0
    SET_TEST_NUM 11
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc12_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 1
    SET_TEST_NUM 12
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 16
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 15
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc14_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, -1
    SET_TEST_NUM 14
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc14_exp, 16
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc15_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, -16
    SET_TEST_NUM 15
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 16
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 0
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc17_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 1
    SET_TEST_NUM 17
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 32
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc18_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, 15
    SET_TEST_NUM 18
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc18_exp, 32
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc19_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, -1
    SET_TEST_NUM 19
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 32
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc20_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vsadd.vi v8, v16, -16
    SET_TEST_NUM 20
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc20_exp, 32
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x00, 0x01, 0xff, 0x7f
tc1_exp:
    .byte 0x00, 0x01, 0xff, 0x7f
.align 1
tc2_s2:
    .byte 0x00, 0x01, 0xff, 0x7f
tc2_exp:
    .byte 0x01, 0x02, 0x00, 0x7f
.align 1
tc3_s2:
    .byte 0x00, 0x01, 0xff, 0x7f
tc3_exp:
    .byte 0x0f, 0x10, 0x0e, 0x7f
.align 1
tc4_s2:
    .byte 0x00, 0x01, 0xff, 0x7f
tc4_exp:
    .byte 0xff, 0x00, 0xfe, 0x7e
.align 1
tc5_s2:
    .byte 0x00, 0x01, 0xff, 0x7f
tc5_exp:
    .byte 0xf0, 0xf1, 0xef, 0x6f
.align 1
tc6_s2:
    .half 0x0000, 0x0001, 0xffff, 0x7fff
tc6_exp:
    .half 0x0000, 0x0001, 0xffff, 0x7fff
.align 1
tc7_s2:
    .half 0x0000, 0x0001, 0xffff, 0x7fff
tc7_exp:
    .half 0x0001, 0x0002, 0x0000, 0x7fff
.align 1
tc8_s2:
    .half 0x0000, 0x0001, 0xffff, 0x7fff
tc8_exp:
    .half 0x000f, 0x0010, 0x000e, 0x7fff
.align 1
tc9_s2:
    .half 0x0000, 0x0001, 0xffff, 0x7fff
tc9_exp:
    .half 0xffff, 0x0000, 0xfffe, 0x7ffe
.align 1
tc10_s2:
    .half 0x0000, 0x0001, 0xffff, 0x7fff
tc10_exp:
    .half 0xfff0, 0xfff1, 0xffef, 0x7fef
.align 2
tc11_s2:
    .word 0x00000000, 0x00000001, 0xffffffff, 0x7fffffff
tc11_exp:
    .word 0x00000000, 0x00000001, 0xffffffff, 0x7fffffff
.align 2
tc12_s2:
    .word 0x00000000, 0x00000001, 0xffffffff, 0x7fffffff
tc12_exp:
    .word 0x00000001, 0x00000002, 0x00000000, 0x7fffffff
.align 2
tc13_s2:
    .word 0x00000000, 0x00000001, 0xffffffff, 0x7fffffff
tc13_exp:
    .word 0x0000000f, 0x00000010, 0x0000000e, 0x7fffffff
.align 2
tc14_s2:
    .word 0x00000000, 0x00000001, 0xffffffff, 0x7fffffff
tc14_exp:
    .word 0xffffffff, 0x00000000, 0xfffffffe, 0x7ffffffe
.align 2
tc15_s2:
    .word 0x00000000, 0x00000001, 0xffffffff, 0x7fffffff
tc15_exp:
    .word 0xfffffff0, 0xfffffff1, 0xffffffef, 0x7fffffef
.align 3
tc16_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0xffffffffffffffff, 0x7fffffffffffffff
tc16_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0xffffffffffffffff, 0x7fffffffffffffff
.align 3
tc17_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0xffffffffffffffff, 0x7fffffffffffffff
tc17_exp:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000000, 0x7fffffffffffffff
.align 3
tc18_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0xffffffffffffffff, 0x7fffffffffffffff
tc18_exp:
    .dword 0x000000000000000f, 0x0000000000000010, 0x000000000000000e, 0x7fffffffffffffff
.align 3
tc19_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0xffffffffffffffff, 0x7fffffffffffffff
tc19_exp:
    .dword 0xffffffffffffffff, 0x0000000000000000, 0xfffffffffffffffe, 0x7ffffffffffffffe
.align 3
tc20_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0xffffffffffffffff, 0x7fffffffffffffff
tc20_exp:
    .dword 0xfffffffffffffff0, 0xfffffffffffffff1, 0xffffffffffffffef, 0x7fffffffffffffef

.align 4
result_buf:  .space 256
witness_buf: .space 256

