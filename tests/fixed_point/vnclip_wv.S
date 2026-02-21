/* Auto-generated test for vnclip.wv
 * Narrowing clip vnclip.wv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vnclip.wv e8 basic: result
 *     2 = vnclip.wv e8 shift: result
 *     3 = vnclip.wv e16 basic: result
 *     4 = vnclip.wv e16 shift: result
 *     5 = vnclip.wv e32 basic: result
 *     6 = vnclip.wv e32 shift: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    csrwi vxrm, 0
    la t1, tc1_s2w
    vle16.v v16, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vnclip.wv v8, v16, v20
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
    la t1, tc2_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vnclip.wv v8, v16, v20
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
    la t1, tc3_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vnclip.wv v8, v16, v20
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
    la t1, tc4_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vnclip.wv v8, v16, v20
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
    la t1, tc5_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vnclip.wv v8, v16, v20
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
    la t1, tc6_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vnclip.wv v8, v16, v20
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
tc1_s1:
    .byte 0x00, 0x00, 0x00, 0x00
tc1_exp:
    .byte 0x7f, 0x7f, 0xff, 0x00
.align 1
tc2_s2w:
    .half 0xffff, 0x7fff, 0x0001, 0x0000
tc2_s1:
    .byte 0x08, 0x04, 0x00, 0x01
tc2_exp:
    .byte 0x00, 0x7f, 0x01, 0x00
.align 2
tc3_s2w:
    .word 0x00000100, 0x00000200, 0xffffffff, 0x00000000
tc3_s1:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc3_exp:
    .half 0x0100, 0x0200, 0xffff, 0x0000
.align 2
tc4_s2w:
    .word 0xffffffff, 0x7fffffff, 0x00000001, 0x00000000
tc4_s1:
    .half 0x0010, 0x0008, 0x0000, 0x0001
tc4_exp:
    .half 0x0000, 0x7fff, 0x0001, 0x0000
.align 3
tc5_s2w:
    .dword 0x0000000000000100, 0x0000000000000200, 0xffffffffffffffff, 0x0000000000000000
tc5_s1:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc5_exp:
    .word 0x00000100, 0x00000200, 0xffffffff, 0x00000000
.align 3
tc6_s2w:
    .dword 0xffffffffffffffff, 0x7fffffffffffffff, 0x0000000000000001, 0x0000000000000000
tc6_s1:
    .word 0x00000020, 0x00000010, 0x00000000, 0x00000001
tc6_exp:
    .word 0x00000000, 0x7fffffff, 0x00000001, 0x00000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

