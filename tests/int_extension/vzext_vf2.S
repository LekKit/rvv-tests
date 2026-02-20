/* Auto-generated test for vzext.vf2
 * Integer extension vzext.vf2
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vzext.vf2 e8→e16 basic: result
 *     2 = vzext.vf2 e8→e16 basic: CSR side-effect
 *     3 = vzext.vf2 e8→e16 max: result
 *     4 = vzext.vf2 e8→e16 max: CSR side-effect
 *     5 = vzext.vf2 e8→e16 high: result
 *     6 = vzext.vf2 e8→e16 high: CSR side-effect
 *     7 = vzext.vf2 e16→e32 basic: result
 *     8 = vzext.vf2 e16→e32 basic: CSR side-effect
 *     9 = vzext.vf2 e16→e32 max: result
 *    10 = vzext.vf2 e16→e32 max: CSR side-effect
 *    11 = vzext.vf2 e16→e32 high: result
 *    12 = vzext.vf2 e16→e32 high: CSR side-effect
 *    13 = vzext.vf2 e32→e64 basic: result
 *    14 = vzext.vf2 e32→e64 basic: CSR side-effect
 *    15 = vzext.vf2 e32→e64 max: result
 *    16 = vzext.vf2 e32→e64 max: CSR side-effect
 *    17 = vzext.vf2 e32→e64 high: result
 *    18 = vzext.vf2 e32→e64 high: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vzext.vf2 e8→e16 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e16, m2, tu, mu
    SAVE_CSRS
    vzext.vf2 v8, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8
    SET_TEST_NUM 2
    CHECK_VSTART_ZERO

    /* Test 3-4: vzext.vf2 e8→e16 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e16, m2, tu, mu
    SAVE_CSRS
    vzext.vf2 v8, v16
    SET_TEST_NUM 3
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    SET_TEST_NUM 4
    CHECK_VSTART_ZERO

    /* Test 5-6: vzext.vf2 e8→e16 high */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e16, m2, tu, mu
    SAVE_CSRS
    vzext.vf2 v8, v16
    SET_TEST_NUM 5
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8
    SET_TEST_NUM 6
    CHECK_VSTART_ZERO

    /* Test 7-8: vzext.vf2 e16→e32 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_src
    vle16.v v16, (t1)
    vsetvli t0, t0, e32, m2, tu, mu
    SAVE_CSRS
    vzext.vf2 v8, v16
    SET_TEST_NUM 7
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 16
    SET_TEST_NUM 8
    CHECK_VSTART_ZERO

    /* Test 9-10: vzext.vf2 e16→e32 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_src
    vle16.v v16, (t1)
    vsetvli t0, t0, e32, m2, tu, mu
    SAVE_CSRS
    vzext.vf2 v8, v16
    SET_TEST_NUM 9
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16
    SET_TEST_NUM 10
    CHECK_VSTART_ZERO

    /* Test 11-12: vzext.vf2 e16→e32 high */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_src
    vle16.v v16, (t1)
    vsetvli t0, t0, e32, m2, tu, mu
    SAVE_CSRS
    vzext.vf2 v8, v16
    SET_TEST_NUM 11
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16
    SET_TEST_NUM 12
    CHECK_VSTART_ZERO

    /* Test 13-14: vzext.vf2 e32→e64 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_src
    vle32.v v16, (t1)
    vsetvli t0, t0, e64, m2, tu, mu
    SAVE_CSRS
    vzext.vf2 v8, v16
    SET_TEST_NUM 13
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 32
    SET_TEST_NUM 14
    CHECK_VSTART_ZERO

    /* Test 15-16: vzext.vf2 e32→e64 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc15_src
    vle32.v v16, (t1)
    vsetvli t0, t0, e64, m2, tu, mu
    SAVE_CSRS
    vzext.vf2 v8, v16
    SET_TEST_NUM 15
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 32
    SET_TEST_NUM 16
    CHECK_VSTART_ZERO

    /* Test 17-18: vzext.vf2 e32→e64 high */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc17_src
    vle32.v v16, (t1)
    vsetvli t0, t0, e64, m2, tu, mu
    SAVE_CSRS
    vzext.vf2 v8, v16
    SET_TEST_NUM 17
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 32
    SET_TEST_NUM 18
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 1
tc1_src:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_exp:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc3_src:
    .byte 0xff, 0xfe, 0x7f, 0x00
tc3_exp:
    .half 0x00ff, 0x00fe, 0x007f, 0x0000
.align 1
tc5_src:
    .byte 0xff, 0x80, 0x81, 0x01
tc5_exp:
    .half 0x00ff, 0x0080, 0x0081, 0x0001
.align 2
tc7_src:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc7_exp:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc9_src:
    .half 0xffff, 0xfffe, 0x7fff, 0x0000
tc9_exp:
    .word 0x0000ffff, 0x0000fffe, 0x00007fff, 0x00000000
.align 2
tc11_src:
    .half 0xffff, 0x8000, 0x8001, 0x0001
tc11_exp:
    .word 0x0000ffff, 0x00008000, 0x00008001, 0x00000001
.align 3
tc13_src:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc13_exp:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc15_src:
    .word 0xffffffff, 0xfffffffe, 0x7fffffff, 0x00000000
tc15_exp:
    .dword 0x00000000ffffffff, 0x00000000fffffffe, 0x000000007fffffff, 0x0000000000000000
.align 3
tc17_src:
    .word 0xffffffff, 0x80000000, 0x80000001, 0x00000001
tc17_exp:
    .dword 0x00000000ffffffff, 0x0000000080000000, 0x0000000080000001, 0x0000000000000001

.align 4
result_buf:  .space 256
witness_buf: .space 256

