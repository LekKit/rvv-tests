/* Auto-generated test for vssrl.vv
 * Scaling shift vssrl.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vssrl.vv e8 shift0: result
 *     2 = vssrl.vv e8 shift1: result
 *     3 = vssrl.vv e16 shift0: result
 *     4 = vssrl.vv e16 shift1: result
 *     5 = vssrl.vv e32 shift0: result
 *     6 = vssrl.vv e32 shift1: result
 *     7 = vssrl.vv e64 shift0: result
 *     8 = vssrl.vv e64 shift1: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vssrl.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc2_s2
    vle8.v v16, (t1)
    la t1, tc2_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vssrl.vv v8, v16, v20
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc3_s2
    vle16.v v16, (t1)
    la t1, tc3_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vssrl.vv v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc4_s2
    vle16.v v16, (t1)
    la t1, tc4_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vssrl.vv v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc5_s2
    vle32.v v16, (t1)
    la t1, tc5_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vssrl.vv v8, v16, v20
    SET_TEST_NUM 5
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc6_s2
    vle32.v v16, (t1)
    la t1, tc6_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vssrl.vv v8, v16, v20
    SET_TEST_NUM 6
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc7_s2
    vle64.v v16, (t1)
    la t1, tc7_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vssrl.vv v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc8_s2
    vle64.v v16, (t1)
    la t1, tc8_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vssrl.vv v8, v16, v20
    SET_TEST_NUM 8
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0xff, 0x01, 0x00, 0x7f
tc1_s1:
    .byte 0x00, 0x00, 0x00, 0x00
tc1_exp:
    .byte 0xff, 0x01, 0x00, 0x7f
.align 1
tc2_s2:
    .byte 0xff, 0x55, 0xaa, 0x01
tc2_s1:
    .byte 0x01, 0x01, 0x01, 0x01
tc2_exp:
    .byte 0x80, 0x2b, 0x55, 0x01
.align 1
tc3_s2:
    .half 0xffff, 0x0001, 0x0000, 0x7fff
tc3_s1:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc3_exp:
    .half 0xffff, 0x0001, 0x0000, 0x7fff
.align 1
tc4_s2:
    .half 0xffff, 0x0055, 0x00aa, 0x0001
tc4_s1:
    .half 0x0001, 0x0001, 0x0001, 0x0001
tc4_exp:
    .half 0x8000, 0x002b, 0x0055, 0x0001
.align 2
tc5_s2:
    .word 0xffffffff, 0x00000001, 0x00000000, 0x7fffffff
tc5_s1:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc5_exp:
    .word 0xffffffff, 0x00000001, 0x00000000, 0x7fffffff
.align 2
tc6_s2:
    .word 0xffffffff, 0x00000055, 0x000000aa, 0x00000001
tc6_s1:
    .word 0x00000001, 0x00000001, 0x00000001, 0x00000001
tc6_exp:
    .word 0x80000000, 0x0000002b, 0x00000055, 0x00000001
.align 3
tc7_s2:
    .dword 0xffffffffffffffff, 0x0000000000000001, 0x0000000000000000, 0x7fffffffffffffff
tc7_s1:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc7_exp:
    .dword 0xffffffffffffffff, 0x0000000000000001, 0x0000000000000000, 0x7fffffffffffffff
.align 3
tc8_s2:
    .dword 0xffffffffffffffff, 0x0000000000000055, 0x00000000000000aa, 0x0000000000000001
tc8_s1:
    .dword 0x0000000000000001, 0x0000000000000001, 0x0000000000000001, 0x0000000000000001
tc8_exp:
    .dword 0x8000000000000000, 0x000000000000002b, 0x0000000000000055, 0x0000000000000001

.align 4
result_buf:  .space 256
witness_buf: .space 256

