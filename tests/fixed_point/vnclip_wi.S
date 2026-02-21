/* Auto-generated test for vnclip.wi
 * Narrowing clip vnclip.wi
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vnclip.wi e8 imm0: result
 *     2 = vnclip.wi e8 imm1: result
 *     3 = vnclip.wi e16 imm0: result
 *     4 = vnclip.wi e16 imm1: result
 *     5 = vnclip.wi e32 imm0: result
 *     6 = vnclip.wi e32 imm1: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    csrwi vxrm, 0
    la t1, tc1_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vnclip.wi v8, v16, 0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    csrwi vxrm, 0
    la t1, tc2_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vnclip.wi v8, v16, 1
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    csrwi vxrm, 0
    la t1, tc3_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vnclip.wi v8, v16, 0
    SET_TEST_NUM 3
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    csrwi vxrm, 0
    la t1, tc4_s2w
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vnclip.wi v8, v16, 1
    SET_TEST_NUM 4
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 8
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    csrwi vxrm, 0
    la t1, tc5_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vnclip.wi v8, v16, 0
    SET_TEST_NUM 5
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 16
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    csrwi vxrm, 0
    la t1, tc6_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vnclip.wi v8, v16, 1
    SET_TEST_NUM 6
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 16
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    PASS_TEST

.data
.align 1
tc1_s2w:
    .half 0x0100, 0x0200, 0xffff, 0x0000
tc1_exp:
    .byte 0x7f, 0x7f, 0xff, 0x00
.align 1
tc2_s2w:
    .half 0xffff, 0x7fff, 0x0001, 0x0000
tc2_exp:
    .byte 0x00, 0x7f, 0x01, 0x00
.align 2
tc3_s2w:
    .word 0x00000100, 0x00000200, 0xffffffff, 0x00000000
tc3_exp:
    .half 0x0100, 0x0200, 0xffff, 0x0000
.align 2
tc4_s2w:
    .word 0xffffffff, 0x7fffffff, 0x00000001, 0x00000000
tc4_exp:
    .half 0x0000, 0x7fff, 0x0001, 0x0000
.align 3
tc5_s2w:
    .dword 0x0000000000000100, 0x0000000000000200, 0xffffffffffffffff, 0x0000000000000000
tc5_exp:
    .word 0x00000100, 0x00000200, 0xffffffff, 0x00000000
.align 3
tc6_s2w:
    .dword 0xffffffffffffffff, 0x7fffffffffffffff, 0x0000000000000001, 0x0000000000000000
tc6_exp:
    .word 0x00000000, 0x7fffffff, 0x00000001, 0x00000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

