/* Auto-generated test for vwsubu.vv
 * Widening vwsubu.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwsubu.vv e8→e16 basic: wrong result
 *     2 = vwsubu.vv e8→e16 basic: CSR side-effect
 *     3 = vwsubu.vv e8→e16 max: wrong result
 *     4 = vwsubu.vv e8→e16 max: CSR side-effect
 *     5 = vwsubu.vv e8→e16 signed: wrong result
 *     6 = vwsubu.vv e8→e16 signed: CSR side-effect
 *     7 = vwsubu.vv e16→e32 basic: wrong result
 *     8 = vwsubu.vv e16→e32 basic: CSR side-effect
 *     9 = vwsubu.vv e16→e32 max: wrong result
 *    10 = vwsubu.vv e16→e32 max: CSR side-effect
 *    11 = vwsubu.vv e16→e32 signed: wrong result
 *    12 = vwsubu.vv e16→e32 signed: CSR side-effect
 *    13 = vwsubu.vv e32→e64 basic: wrong result
 *    14 = vwsubu.vv e32→e64 basic: CSR side-effect
 *    15 = vwsubu.vv e32→e64 max: wrong result
 *    16 = vwsubu.vv e32→e64 max: CSR side-effect
 *    17 = vwsubu.vv e32→e64 signed: wrong result
 *    18 = vwsubu.vv e32→e64 signed: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vwsubu.vv SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwsubu.vv v8, v16, v20
    SET_TEST_NUM 1
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8
    SET_TEST_NUM 2
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_FCSR_UNCHANGED

    /* Test 3-4: vwsubu.vv SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    la t1, tc3_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwsubu.vv v8, v16, v20
    SET_TEST_NUM 3
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    SET_TEST_NUM 4
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_FCSR_UNCHANGED

    /* Test 5-6: vwsubu.vv SEW=8 signed */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_s2
    vle8.v v16, (t1)
    la t1, tc5_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vwsubu.vv v8, v16, v20
    SET_TEST_NUM 5
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8
    SET_TEST_NUM 6
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_FCSR_UNCHANGED

    /* Test 7-8: vwsubu.vv SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_s2
    vle16.v v16, (t1)
    la t1, tc7_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwsubu.vv v8, v16, v20
    SET_TEST_NUM 7
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 16
    SET_TEST_NUM 8
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_FCSR_UNCHANGED

    /* Test 9-10: vwsubu.vv SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_s2
    vle16.v v16, (t1)
    la t1, tc9_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwsubu.vv v8, v16, v20
    SET_TEST_NUM 9
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16
    SET_TEST_NUM 10
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_FCSR_UNCHANGED

    /* Test 11-12: vwsubu.vv SEW=16 signed */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_s2
    vle16.v v16, (t1)
    la t1, tc11_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vwsubu.vv v8, v16, v20
    SET_TEST_NUM 11
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16
    SET_TEST_NUM 12
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_FCSR_UNCHANGED

    /* Test 13-14: vwsubu.vv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_s2
    vle32.v v16, (t1)
    la t1, tc13_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwsubu.vv v8, v16, v20
    SET_TEST_NUM 13
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 32
    SET_TEST_NUM 14
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_FCSR_UNCHANGED

    /* Test 15-16: vwsubu.vv SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc15_s2
    vle32.v v16, (t1)
    la t1, tc15_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwsubu.vv v8, v16, v20
    SET_TEST_NUM 15
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 32
    SET_TEST_NUM 16
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_FCSR_UNCHANGED

    /* Test 17-18: vwsubu.vv SEW=32 signed */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc17_s2
    vle32.v v16, (t1)
    la t1, tc17_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwsubu.vv v8, v16, v20
    SET_TEST_NUM 17
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 32
    SET_TEST_NUM 18
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_FCSR_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_s1:
    .byte 0x05, 0x06, 0x07, 0x08
tc1_exp:
    .half 0xfffc, 0xfffc, 0xfffc, 0xfffc
.align 1
tc3_s2:
    .byte 0xff, 0xff, 0xff, 0xff
tc3_s1:
    .byte 0xff, 0xff, 0xff, 0xff
tc3_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 1
tc5_s2:
    .byte 0x80, 0x81, 0x00, 0xff
tc5_s1:
    .byte 0x01, 0xff, 0x80, 0x00
tc5_exp:
    .half 0x007f, 0xff82, 0xff80, 0x00ff
.align 2
tc7_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc7_s1:
    .half 0x0005, 0x0006, 0x0007, 0x0008
tc7_exp:
    .word 0xfffffffc, 0xfffffffc, 0xfffffffc, 0xfffffffc
.align 2
tc9_s2:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc9_s1:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc9_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
.align 2
tc11_s2:
    .half 0x8000, 0x8001, 0x0000, 0xffff
tc11_s1:
    .half 0x0001, 0xffff, 0x8000, 0x0000
tc11_exp:
    .word 0x00007fff, 0xffff8002, 0xffff8000, 0x0000ffff
.align 3
tc13_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc13_s1:
    .word 0x00000005, 0x00000006, 0x00000007, 0x00000008
tc13_exp:
    .dword 0xfffffffffffffffc, 0xfffffffffffffffc, 0xfffffffffffffffc, 0xfffffffffffffffc
.align 3
tc15_s2:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc15_s1:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc15_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc17_s2:
    .word 0x80000000, 0x80000001, 0x00000000, 0xffffffff
tc17_s1:
    .word 0x00000001, 0xffffffff, 0x80000000, 0x00000000
tc17_exp:
    .dword 0x000000007fffffff, 0xffffffff80000002, 0xffffffff80000000, 0x00000000ffffffff

.align 4
result_buf:  .space 256
witness_buf: .space 256

