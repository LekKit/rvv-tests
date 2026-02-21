/* Auto-generated test for vssra.vx
 * Scaling shift vssra.vx
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vssra.vx e8 sh=0: result
 *     2 = vssra.vx e8 sh=1: result
 *     3 = vssra.vx e8 sh=7: result
 *     4 = vssra.vx e16 sh=0: result
 *     5 = vssra.vx e16 sh=1: result
 *     6 = vssra.vx e16 sh=15: result
 *     7 = vssra.vx e32 sh=0: result
 *     8 = vssra.vx e32 sh=1: result
 *     9 = vssra.vx e32 sh=31: result
 *    10 = vssra.vx e64 sh=0: result
 *    11 = vssra.vx e64 sh=1: result
 *    12 = vssra.vx e64 sh=63: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc1_s2
    vle8.v v16, (t1)
    li a0, 0
    SAVE_CSRS
    vssra.vx v8, v16, a0
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
    li a0, 1
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc3_s2
    vle8.v v16, (t1)
    li a0, 7
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc4_s2
    vle16.v v16, (t1)
    li a0, 0
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 4
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc5_s2
    vle16.v v16, (t1)
    li a0, 1
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 5
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc6_s2
    vle16.v v16, (t1)
    li a0, 15
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 6
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc7_s2
    vle32.v v16, (t1)
    li a0, 0
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 7
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc8_s2
    vle32.v v16, (t1)
    li a0, 1
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 8
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc9_s2
    vle32.v v16, (t1)
    li a0, 31
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 9
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc10_s2
    vle64.v v16, (t1)
    li a0, 0
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 10
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc11_s2
    vle64.v v16, (t1)
    li a0, 1
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 11
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc12_s2
    vle64.v v16, (t1)
    li a0, 63
    SAVE_CSRS
    vssra.vx v8, v16, a0
    SET_TEST_NUM 12
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0xff, 0x55, 0xaa, 0x01
tc1_exp:
    .byte 0xff, 0x55, 0xaa, 0x01
.align 1
tc2_s2:
    .byte 0xff, 0x55, 0xaa, 0x01
tc2_exp:
    .byte 0x00, 0x2b, 0xd5, 0x01
.align 1
tc3_s2:
    .byte 0xff, 0x55, 0xaa, 0x01
tc3_exp:
    .byte 0x00, 0x01, 0xff, 0x00
.align 1
tc4_s2:
    .half 0xffff, 0x0055, 0x00aa, 0x0001
tc4_exp:
    .half 0xffff, 0x0055, 0x00aa, 0x0001
.align 1
tc5_s2:
    .half 0xffff, 0x0055, 0x00aa, 0x0001
tc5_exp:
    .half 0x0000, 0x002b, 0x0055, 0x0001
.align 1
tc6_s2:
    .half 0xffff, 0x0055, 0x00aa, 0x0001
tc6_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 2
tc7_s2:
    .word 0xffffffff, 0x00000055, 0x000000aa, 0x00000001
tc7_exp:
    .word 0xffffffff, 0x00000055, 0x000000aa, 0x00000001
.align 2
tc8_s2:
    .word 0xffffffff, 0x00000055, 0x000000aa, 0x00000001
tc8_exp:
    .word 0x00000000, 0x0000002b, 0x00000055, 0x00000001
.align 2
tc9_s2:
    .word 0xffffffff, 0x00000055, 0x000000aa, 0x00000001
tc9_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
.align 3
tc10_s2:
    .dword 0xffffffffffffffff, 0x0000000000000055, 0x00000000000000aa, 0x0000000000000001
tc10_exp:
    .dword 0xffffffffffffffff, 0x0000000000000055, 0x00000000000000aa, 0x0000000000000001
.align 3
tc11_s2:
    .dword 0xffffffffffffffff, 0x0000000000000055, 0x00000000000000aa, 0x0000000000000001
tc11_exp:
    .dword 0x0000000000000000, 0x000000000000002b, 0x0000000000000055, 0x0000000000000001
.align 3
tc12_s2:
    .dword 0xffffffffffffffff, 0x0000000000000055, 0x00000000000000aa, 0x0000000000000001
tc12_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

