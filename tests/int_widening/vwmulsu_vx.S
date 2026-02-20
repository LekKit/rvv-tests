/* Auto-generated test for vwmulsu.vx
 * Widening vwmulsu.vx
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwmulsu.vx e8 sc_0x01: result
 *     2 = vwmulsu.vx e8 sc_0xff: result
 *     3 = vwmulsu.vx e8 sc_0x00: result
 *     4 = vwmulsu.vx e8 sc_0x80: result
 *     5 = vwmulsu.vx e16 sc_0x0001: result
 *     6 = vwmulsu.vx e16 sc_0xffff: result
 *     7 = vwmulsu.vx e16 sc_0x0000: result
 *     8 = vwmulsu.vx e16 sc_0x8000: result
 *     9 = vwmulsu.vx e32 sc_0x00000001: result
 *    10 = vwmulsu.vx e32 sc_0xffffffff: result
 *    11 = vwmulsu.vx e32 sc_0x00000000: result
 *    12 = vwmulsu.vx e32 sc_0x80000000: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vwmulsu.vx SEW=8 sc_0x01 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    li a0, 0x01
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 1
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8
    CHECK_VSTART_ZERO

    /* Test 2: vwmulsu.vx SEW=8 sc_0xff */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_s2
    vle8.v v16, (t1)
    li a0, 0xff
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 2
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_VSTART_ZERO

    /* Test 3: vwmulsu.vx SEW=8 sc_0x00 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    li a0, 0x00
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 3
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    CHECK_VSTART_ZERO

    /* Test 4: vwmulsu.vx SEW=8 sc_0x80 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    li a0, 0x80
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 8
    CHECK_VSTART_ZERO

    /* Test 5: vwmulsu.vx SEW=16 sc_0x0001 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_s2
    vle16.v v16, (t1)
    li a0, 0x0001
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 5
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 16
    CHECK_VSTART_ZERO

    /* Test 6: vwmulsu.vx SEW=16 sc_0xffff */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_s2
    vle16.v v16, (t1)
    li a0, 0xffff
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 6
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 16
    CHECK_VSTART_ZERO

    /* Test 7: vwmulsu.vx SEW=16 sc_0x0000 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_s2
    vle16.v v16, (t1)
    li a0, 0x0000
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 7
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 16
    CHECK_VSTART_ZERO

    /* Test 8: vwmulsu.vx SEW=16 sc_0x8000 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc8_s2
    vle16.v v16, (t1)
    li a0, 0x8000
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 8
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 16
    CHECK_VSTART_ZERO

    /* Test 9: vwmulsu.vx SEW=32 sc_0x00000001 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_s2
    vle32.v v16, (t1)
    li a0, 0x00000001
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 9
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 32
    CHECK_VSTART_ZERO

    /* Test 10: vwmulsu.vx SEW=32 sc_0xffffffff */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc10_s2
    vle32.v v16, (t1)
    li a0, 0xffffffff
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 10
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 32
    CHECK_VSTART_ZERO

    /* Test 11: vwmulsu.vx SEW=32 sc_0x00000000 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_s2
    vle32.v v16, (t1)
    li a0, 0x00000000
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 11
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 32
    CHECK_VSTART_ZERO

    /* Test 12: vwmulsu.vx SEW=32 sc_0x80000000 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc12_s2
    vle32.v v16, (t1)
    li a0, 0x80000000
    SAVE_CSRS
    vwmulsu.vx v8, v16, a0
    SET_TEST_NUM 12
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 32
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0xff, 0x7f
tc1_exp:
    .half 0x0001, 0x0002, 0xffff, 0x007f
.align 1
tc2_s2:
    .byte 0x01, 0x02, 0xff, 0x7f
tc2_exp:
    .half 0x00ff, 0x01fe, 0xff01, 0x7e81
.align 1
tc3_s2:
    .byte 0x01, 0x02, 0xff, 0x7f
tc3_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 1
tc4_s2:
    .byte 0x01, 0x02, 0xff, 0x7f
tc4_exp:
    .half 0x0080, 0x0100, 0xff80, 0x3f80
.align 2
tc5_s2:
    .half 0x0001, 0x0002, 0xffff, 0x7fff
tc5_exp:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x00007fff
.align 2
tc6_s2:
    .half 0x0001, 0x0002, 0xffff, 0x7fff
tc6_exp:
    .word 0x0000ffff, 0x0001fffe, 0xffff0001, 0x7ffe8001
.align 2
tc7_s2:
    .half 0x0001, 0x0002, 0xffff, 0x7fff
tc7_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
.align 2
tc8_s2:
    .half 0x0001, 0x0002, 0xffff, 0x7fff
tc8_exp:
    .word 0x00008000, 0x00010000, 0xffff8000, 0x3fff8000
.align 3
tc9_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x7fffffff
tc9_exp:
    .dword 0x0000000000000001, 0x0000000000000002, 0xffffffffffffffff, 0x000000007fffffff
.align 3
tc10_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x7fffffff
tc10_exp:
    .dword 0x00000000ffffffff, 0x00000001fffffffe, 0xffffffff00000001, 0x7ffffffe80000001
.align 3
tc11_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x7fffffff
tc11_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc12_s2:
    .word 0x00000001, 0x00000002, 0xffffffff, 0x7fffffff
tc12_exp:
    .dword 0x0000000080000000, 0x0000000100000000, 0xffffffff80000000, 0x3fffffff80000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

