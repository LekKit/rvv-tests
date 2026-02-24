/* Auto-generated test for vsext.vf4
 * Integer extension vsext.vf4
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsext.vf4 e8→e32 basic: result
 *     2 = vsext.vf4 e8→e32 basic: CSR side-effect
 *     3 = vsext.vf4 e8→e32 max: result
 *     4 = vsext.vf4 e8→e32 max: CSR side-effect
 *     5 = vsext.vf4 e8→e32 high: result
 *     6 = vsext.vf4 e8→e32 high: CSR side-effect
 *     7 = vsext.vf4 e16→e64 basic: result
 *     8 = vsext.vf4 e16→e64 basic: CSR side-effect
 *     9 = vsext.vf4 e16→e64 max: result
 *    10 = vsext.vf4 e16→e64 max: CSR side-effect
 *    11 = vsext.vf4 e16→e64 high: result
 *    12 = vsext.vf4 e16→e64 high: CSR side-effect
 *    13 = vsext.vf4 e8→e32 masked: result
 *    14 = vsext.vf4 e8→e32 masked: witness
 *    15 = vsext.vf4 e8→e32 masked: CSR side-effect
 *    16 = vsext.vf4 e16→e64 masked: result
 *    17 = vsext.vf4 e16→e64 masked: witness
 *    18 = vsext.vf4 e16→e64 masked: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vsext.vf4 e8→e32 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e32, m4, tu, mu
    SAVE_CSRS
    vsext.vf4 v8, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    /* Test 3-4: vsext.vf4 e8→e32 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e32, m4, tu, mu
    SAVE_CSRS
    vsext.vf4 v8, v16
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 16
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED

    /* Test 5-6: vsext.vf4 e8→e32 high */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e32, m4, tu, mu
    SAVE_CSRS
    vsext.vf4 v8, v16
    SET_TEST_NUM 5
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 16
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-8: vsext.vf4 e16→e64 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_src
    vle16.v v16, (t1)
    vsetvli t0, t0, e64, m4, tu, mu
    SAVE_CSRS
    vsext.vf4 v8, v16
    SET_TEST_NUM 7
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 32
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED

    /* Test 9-10: vsext.vf4 e16→e64 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_src
    vle16.v v16, (t1)
    vsetvli t0, t0, e64, m4, tu, mu
    SAVE_CSRS
    vsext.vf4 v8, v16
    SET_TEST_NUM 9
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 32
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED

    /* Test 11-12: vsext.vf4 e16→e64 high */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_src
    vle16.v v16, (t1)
    vsetvli t0, t0, e64, m4, tu, mu
    SAVE_CSRS
    vsext.vf4 v8, v16
    SET_TEST_NUM 11
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 32
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13-15: vsext.vf4 e8→e32 masked (v0.t, mu) */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc13_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e32, m4, tu, mu
    la t1, tc13_vd
    vle32.v v8, (t1)
    la t1, tc13_mask
    vlm.v v0, (t1)
    la t1, tc13_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vsext.vf4 v8, v16, v0.t
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16
    SET_TEST_NUM 14
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 16
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED

    /* Test 16-18: vsext.vf4 e16→e64 masked (v0.t, mu) */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc16_src
    vle16.v v16, (t1)
    vsetvli t0, t0, e64, m4, tu, mu
    la t1, tc16_vd
    vle64.v v8, (t1)
    la t1, tc16_mask
    vlm.v v0, (t1)
    la t1, tc16_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vsext.vf4 v8, v16, v0.t
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    SET_TEST_NUM 17
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc16_w, 32
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_src:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_exp:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc3_src:
    .byte 0xff, 0xfe, 0x7f, 0x00
tc3_exp:
    .word 0xffffffff, 0xfffffffe, 0x0000007f, 0x00000000
.align 2
tc5_src:
    .byte 0xff, 0x80, 0x81, 0x01
tc5_exp:
    .word 0xffffffff, 0xffffff80, 0xffffff81, 0x00000001
.align 3
tc7_src:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc7_exp:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc9_src:
    .half 0xffff, 0xfffe, 0x7fff, 0x0000
tc9_exp:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000007fff, 0x0000000000000000
.align 3
tc11_src:
    .half 0xffff, 0x8000, 0x8001, 0x0001
tc11_exp:
    .dword 0xffffffffffffffff, 0xffffffffffff8000, 0xffffffffffff8001, 0x0000000000000001
.align 1
tc13_mask:
    .byte 5
.align 2
tc13_src:
    .byte 0xff, 0x7f, 0x80, 0x01
tc13_vd:
    .word 0x000000dd, 0x000000dd, 0x000000dd, 0x000000dd
tc13_exp:
    .word 0xffffffff, 0x000000dd, 0xffffff80, 0x000000dd
tc13_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 1
tc16_mask:
    .byte 5
.align 3
tc16_src:
    .half 0xffff, 0x7fff, 0x8000, 0x0001
tc16_vd:
    .dword 0x00000000000000dd, 0x00000000000000dd, 0x00000000000000dd, 0x00000000000000dd
tc16_exp:
    .dword 0xffffffffffffffff, 0x00000000000000dd, 0xffffffffffff8000, 0x00000000000000dd
tc16_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

