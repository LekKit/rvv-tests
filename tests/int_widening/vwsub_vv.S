/* Auto-generated test for vwsub.vv
 * Widening vwsub.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwsub.vv e8→e16 basic: wrong result
 *     2 = vwsub.vv e8→e16 basic: witness changed
 *     3 = vwsub.vv e8→e16 basic: CSR side-effect
 *     4 = vwsub.vv e8→e16 max: wrong result
 *     5 = vwsub.vv e8→e16 max: witness changed
 *     6 = vwsub.vv e8→e16 max: CSR side-effect
 *     7 = vwsub.vv e8→e16 signed: wrong result
 *     8 = vwsub.vv e8→e16 signed: witness changed
 *     9 = vwsub.vv e8→e16 signed: CSR side-effect
 *    10 = vwsub.vv e16→e32 basic: wrong result
 *    11 = vwsub.vv e16→e32 basic: witness changed
 *    12 = vwsub.vv e16→e32 basic: CSR side-effect
 *    13 = vwsub.vv e16→e32 max: wrong result
 *    14 = vwsub.vv e16→e32 max: witness changed
 *    15 = vwsub.vv e16→e32 max: CSR side-effect
 *    16 = vwsub.vv e16→e32 signed: wrong result
 *    17 = vwsub.vv e16→e32 signed: witness changed
 *    18 = vwsub.vv e16→e32 signed: CSR side-effect
 *    19 = vwsub.vv e16→e32 basic masked: wrong result
 *    20 = vwsub.vv e16→e32 basic masked: witness changed
 *    21 = vwsub.vv e16→e32 basic masked: CSR side-effect
 *    22 = vwsub.vv e32→e64 basic: wrong result
 *    23 = vwsub.vv e32→e64 basic: witness changed
 *    24 = vwsub.vv e32→e64 basic: CSR side-effect
 *    25 = vwsub.vv e32→e64 max: wrong result
 *    26 = vwsub.vv e32→e64 max: witness changed
 *    27 = vwsub.vv e32→e64 max: CSR side-effect
 *    28 = vwsub.vv e32→e64 signed: wrong result
 *    29 = vwsub.vv e32→e64 signed: witness changed
 *    30 = vwsub.vv e32→e64 signed: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vwsub.vv SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    la t1, tc1_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vwsub.vv v8, v16, v20
    SET_TEST_NUM 2
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 4
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8

    /* Test 4-6: vwsub.vv SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    la t1, tc4_s1
    vle8.v v20, (t1)
    la t1, tc4_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vwsub.vv v8, v16, v20
    SET_TEST_NUM 5
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 4
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 4
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 8

    /* Test 7-9: vwsub.vv SEW=8 signed */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_s2
    vle8.v v16, (t1)
    la t1, tc7_s1
    vle8.v v20, (t1)
    la t1, tc7_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vwsub.vv v8, v16, v20
    SET_TEST_NUM 8
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 4
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 7
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8

    /* Test 10-12: vwsub.vv SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc10_s2
    vle16.v v16, (t1)
    la t1, tc10_s1
    vle16.v v20, (t1)
    la t1, tc10_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vwsub.vv v8, v16, v20
    SET_TEST_NUM 11
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 8
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 10
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 16

    /* Test 13-15: vwsub.vv SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc13_s2
    vle16.v v16, (t1)
    la t1, tc13_s1
    vle16.v v20, (t1)
    la t1, tc13_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vwsub.vv v8, v16, v20
    SET_TEST_NUM 14
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 8
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 13
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16

    /* Test 16-18: vwsub.vv SEW=16 signed */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc16_s2
    vle16.v v16, (t1)
    la t1, tc16_s1
    vle16.v v20, (t1)
    la t1, tc16_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vwsub.vv v8, v16, v20
    SET_TEST_NUM 17
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc16_w, 8
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 16
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 16

    /* Test 19-21: vwsub.vv SEW=16 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc19_vd
    vle32.v v8, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc19_s2
    vle16.v v16, (t1)
    la t1, tc19_s1
    vle16.v v20, (t1)
    la t1, tc19_mask
    vlm.v v0, (t1)
    la t1, tc19_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vwsub.vv v8, v16, v20, v0.t
    SET_TEST_NUM 20
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc19_w, 8
    SET_TEST_NUM 21
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 19
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 16

    /* Test 22-24: vwsub.vv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc22_s2
    vle32.v v16, (t1)
    la t1, tc22_s1
    vle32.v v20, (t1)
    la t1, tc22_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vwsub.vv v8, v16, v20
    SET_TEST_NUM 23
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc22_w, 16
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 22
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc22_exp, 32

    /* Test 25-27: vwsub.vv SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc25_s2
    vle32.v v16, (t1)
    la t1, tc25_s1
    vle32.v v20, (t1)
    la t1, tc25_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vwsub.vv v8, v16, v20
    SET_TEST_NUM 26
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc25_w, 16
    SET_TEST_NUM 27
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 25
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc25_exp, 32

    /* Test 28-30: vwsub.vv SEW=32 signed */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc28_s2
    vle32.v v16, (t1)
    la t1, tc28_s1
    vle32.v v20, (t1)
    la t1, tc28_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vwsub.vv v8, v16, v20
    SET_TEST_NUM 29
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc28_w, 16
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 28
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc28_exp, 32

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_s1:
    .byte 0x05, 0x06, 0x07, 0x08
tc1_exp:
    .half 0xfffc, 0xfffc, 0xfffc, 0xfffc
tc1_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc4_s2:
    .byte 0xff, 0xff, 0xff, 0xff
tc4_s1:
    .byte 0xff, 0xff, 0xff, 0xff
tc4_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc4_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc7_s2:
    .byte 0x80, 0x81, 0x00, 0xff
tc7_s1:
    .byte 0x01, 0xff, 0x80, 0x00
tc7_exp:
    .half 0xff7f, 0xff82, 0x0080, 0xffff
tc7_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 2
tc10_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc10_s1:
    .half 0x0005, 0x0006, 0x0007, 0x0008
tc10_exp:
    .word 0xfffffffc, 0xfffffffc, 0xfffffffc, 0xfffffffc
tc10_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 2
tc13_s2:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc13_s1:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc13_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc13_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 2
tc16_s2:
    .half 0x8000, 0x8001, 0x0000, 0xffff
tc16_s1:
    .half 0x0001, 0xffff, 0x8000, 0x0000
tc16_exp:
    .word 0xffff7fff, 0xffff8002, 0x00008000, 0xffffffff
tc16_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc19_mask:
    .byte 10
.align 2
tc19_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc19_s1:
    .half 0x0005, 0x0006, 0x0007, 0x0008
tc19_vd:
    .word 0x0000dead, 0x0000dead, 0x0000dead, 0x0000dead
tc19_exp:
    .word 0x0000dead, 0xfffffffc, 0x0000dead, 0xfffffffc
tc19_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 3
tc22_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc22_s1:
    .word 0x00000005, 0x00000006, 0x00000007, 0x00000008
tc22_exp:
    .dword 0xfffffffffffffffc, 0xfffffffffffffffc, 0xfffffffffffffffc, 0xfffffffffffffffc
tc22_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc25_s2:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc25_s1:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc25_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc25_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc28_s2:
    .word 0x80000000, 0x80000001, 0x00000000, 0xffffffff
tc28_s1:
    .word 0x00000001, 0xffffffff, 0x80000000, 0x00000000
tc28_exp:
    .dword 0xffffffff7fffffff, 0xffffffff80000002, 0x0000000080000000, 0xffffffffffffffff
tc28_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0

.align 4
result_buf:  .space 256
witness_buf: .space 256

