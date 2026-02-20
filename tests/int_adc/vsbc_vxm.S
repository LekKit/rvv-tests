/* Auto-generated test for vsbc.vxm
 * Subtract with borrow (vxm)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsbc.vxm e8 basic: result
 *     2 = vsbc.vxm e8 carry: result
 *     3 = vsbc.vxm e8 overflow: result
 *     4 = vsbc.vxm e8 mixed: result
 *     5 = vsbc.vxm e16 basic: result
 *     6 = vsbc.vxm e16 carry: result
 *     7 = vsbc.vxm e16 overflow: result
 *     8 = vsbc.vxm e16 mixed: result
 *     9 = vsbc.vxm e32 basic: result
 *    10 = vsbc.vxm e32 carry: result
 *    11 = vsbc.vxm e32 overflow: result
 *    12 = vsbc.vxm e32 mixed: result
 *    13 = vsbc.vxm e64 basic: result
 *    14 = vsbc.vxm e64 carry: result
 *    15 = vsbc.vxm e64 overflow: result
 *    16 = vsbc.vxm e64 mixed: result
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    vsbc.vxm v8, v16, a0, v0
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
    .byte 0xf7, 0xf8, 0xf9, 0xfa
.align 1
tc2_mask:
    .byte 15
.align 1
tc2_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc2_exp:
    .byte 0xf6, 0xf7, 0xf8, 0xf9
.align 1
tc3_mask:
    .byte 10
.align 1
tc3_s2:
    .byte 0xff, 0xff, 0xfe, 0x00
tc3_exp:
    .byte 0xfe, 0xfd, 0xfd, 0xfe
.align 1
tc4_mask:
    .byte 5
.align 1
tc4_s2:
    .byte 0x00, 0xff, 0x01, 0x7f
tc4_exp:
    .byte 0x00, 0x00, 0x01, 0x80
.align 1
tc5_mask:
    .byte 0
.align 1
tc5_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc5_exp:
    .half 0xfff7, 0xfff8, 0xfff9, 0xfffa
.align 1
tc6_mask:
    .byte 15
.align 1
tc6_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc6_exp:
    .half 0xfff6, 0xfff7, 0xfff8, 0xfff9
.align 1
tc7_mask:
    .byte 10
.align 1
tc7_s2:
    .half 0xffff, 0xffff, 0xfffe, 0x0000
tc7_exp:
    .half 0xfffe, 0xfffd, 0xfffd, 0xfffe
.align 1
tc8_mask:
    .byte 5
.align 1
tc8_s2:
    .half 0x0000, 0xffff, 0x0001, 0x7fff
tc8_exp:
    .half 0x0000, 0x0000, 0x0001, 0x8000
.align 1
tc9_mask:
    .byte 0
.align 2
tc9_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc9_exp:
    .word 0xfffffff7, 0xfffffff8, 0xfffffff9, 0xfffffffa
.align 1
tc10_mask:
    .byte 15
.align 2
tc10_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc10_exp:
    .word 0xfffffff6, 0xfffffff7, 0xfffffff8, 0xfffffff9
.align 1
tc11_mask:
    .byte 10
.align 2
tc11_s2:
    .word 0xffffffff, 0xffffffff, 0xfffffffe, 0x00000000
tc11_exp:
    .word 0xfffffffe, 0xfffffffd, 0xfffffffd, 0xfffffffe
.align 1
tc12_mask:
    .byte 5
.align 2
tc12_s2:
    .word 0x00000000, 0xffffffff, 0x00000001, 0x7fffffff
tc12_exp:
    .word 0x00000000, 0x00000000, 0x00000001, 0x80000000
.align 1
tc13_mask:
    .byte 0
.align 3
tc13_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc13_exp:
    .dword 0xfffffffffffffff7, 0xfffffffffffffff8, 0xfffffffffffffff9, 0xfffffffffffffffa
.align 1
tc14_mask:
    .byte 15
.align 3
tc14_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc14_exp:
    .dword 0xfffffffffffffff6, 0xfffffffffffffff7, 0xfffffffffffffff8, 0xfffffffffffffff9
.align 1
tc15_mask:
    .byte 10
.align 3
tc15_s2:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000
tc15_exp:
    .dword 0xfffffffffffffffe, 0xfffffffffffffffd, 0xfffffffffffffffd, 0xfffffffffffffffe
.align 1
tc16_mask:
    .byte 5
.align 3
tc16_s2:
    .dword 0x0000000000000000, 0xffffffffffffffff, 0x0000000000000001, 0x7fffffffffffffff
tc16_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000001, 0x8000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

