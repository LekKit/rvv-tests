/* Auto-generated test for vaadd.vx
 * Averaging vaadd.vx
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vaadd.vx e8 basic: result
 *     2 = vaadd.vx e8 zero_sc: result
 *     3 = vaadd.vx e8 max_sc: result
 *     4 = vaadd.vx e8 max_vec: result
 *     5 = vaadd.vx e8 signed_sc: result
 *     6 = vaadd.vx e16 basic: result
 *     7 = vaadd.vx e16 zero_sc: result
 *     8 = vaadd.vx e16 max_sc: result
 *     9 = vaadd.vx e16 max_vec: result
 *    10 = vaadd.vx e16 signed_sc: result
 *    11 = vaadd.vx e32 basic: result
 *    12 = vaadd.vx e32 zero_sc: result
 *    13 = vaadd.vx e32 max_sc: result
 *    14 = vaadd.vx e32 max_vec: result
 *    15 = vaadd.vx e32 signed_sc: result
 *    16 = vaadd.vx e64 basic: result
 *    17 = vaadd.vx e64 zero_sc: result
 *    18 = vaadd.vx e64 max_sc: result
 *    19 = vaadd.vx e64 max_vec: result
 *    20 = vaadd.vx e64 signed_sc: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc1_s2
    vle8.v v16, (t1)
    li a0, 0x0a
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc2_s2
    vle8.v v16, (t1)
    li a0, 0x00
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc3_s2
    vle8.v v16, (t1)
    li a0, 0xff
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc4_s2
    vle8.v v16, (t1)
    li a0, 0x01
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc5_s2
    vle8.v v16, (t1)
    li a0, 0x80
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 5
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 4
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc6_s2
    vle16.v v16, (t1)
    li a0, 0x000a
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 6
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 8
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc7_s2
    vle16.v v16, (t1)
    li a0, 0x0000
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 7
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc8_s2
    vle16.v v16, (t1)
    li a0, 0xffff
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 8
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 8
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc9_s2
    vle16.v v16, (t1)
    li a0, 0x0001
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 9
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 8
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc10_s2
    vle16.v v16, (t1)
    li a0, 0x8000
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 10
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 8
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc11_s2
    vle32.v v16, (t1)
    li a0, 0x0000000a
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 11
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc12_s2
    vle32.v v16, (t1)
    li a0, 0x00000000
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 12
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 16
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc13_s2
    vle32.v v16, (t1)
    li a0, 0xffffffff
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc14_s2
    vle32.v v16, (t1)
    li a0, 0x00000001
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 14
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc14_exp, 16
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc15_s2
    vle32.v v16, (t1)
    li a0, 0x80000000
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 15
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 16
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc16_s2
    vle64.v v16, (t1)
    li a0, 0x000000000000000a
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc17_s2
    vle64.v v16, (t1)
    li a0, 0x0000000000000000
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 17
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 32
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc18_s2
    vle64.v v16, (t1)
    li a0, 0xffffffffffffffff
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 18
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc18_exp, 32
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc19_s2
    vle64.v v16, (t1)
    li a0, 0x0000000000000001
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 19
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 32
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    csrwi vxrm, 0
    la t1, tc20_s2
    vle64.v v16, (t1)
    li a0, 0x8000000000000000
    SAVE_CSRS
    vaadd.vx v8, v16, a0
    SET_TEST_NUM 20
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc20_exp, 32
    CHECK_VSTART_ZERO
    CHECK_VXRM_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_exp:
    .byte 0x06, 0x06, 0x07, 0x07
.align 1
tc2_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc2_exp:
    .byte 0x01, 0x01, 0x02, 0x02
.align 1
tc3_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc3_exp:
    .byte 0x00, 0x01, 0x01, 0x02
.align 1
tc4_s2:
    .byte 0xff, 0xff, 0xff, 0xff
tc4_exp:
    .byte 0x00, 0x00, 0x00, 0x00
.align 1
tc5_s2:
    .byte 0x00, 0x01, 0x80, 0xff
tc5_exp:
    .byte 0xc0, 0xc1, 0x80, 0xc0
.align 1
tc6_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc6_exp:
    .half 0x0006, 0x0006, 0x0007, 0x0007
.align 1
tc7_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc7_exp:
    .half 0x0001, 0x0001, 0x0002, 0x0002
.align 1
tc8_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc8_exp:
    .half 0x0000, 0x0001, 0x0001, 0x0002
.align 1
tc9_s2:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc9_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 1
tc10_s2:
    .half 0x0000, 0x0001, 0x8000, 0xffff
tc10_exp:
    .half 0xc000, 0xc001, 0x8000, 0xc000
.align 2
tc11_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc11_exp:
    .word 0x00000006, 0x00000006, 0x00000007, 0x00000007
.align 2
tc12_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc12_exp:
    .word 0x00000001, 0x00000001, 0x00000002, 0x00000002
.align 2
tc13_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc13_exp:
    .word 0x00000000, 0x00000001, 0x00000001, 0x00000002
.align 2
tc14_s2:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc14_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
.align 2
tc15_s2:
    .word 0x00000000, 0x00000001, 0x80000000, 0xffffffff
tc15_exp:
    .word 0xc0000000, 0xc0000001, 0x80000000, 0xc0000000
.align 3
tc16_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc16_exp:
    .dword 0x0000000000000006, 0x0000000000000006, 0x0000000000000007, 0x0000000000000007
.align 3
tc17_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc17_exp:
    .dword 0x0000000000000001, 0x0000000000000001, 0x0000000000000002, 0x0000000000000002
.align 3
tc18_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc18_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0x0000000000000001, 0x0000000000000002
.align 3
tc19_s2:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
tc19_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc20_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0x8000000000000000, 0xffffffffffffffff
tc20_exp:
    .dword 0xc000000000000000, 0xc000000000000001, 0x8000000000000000, 0xc000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

