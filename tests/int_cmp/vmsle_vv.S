/* Auto-generated test for vmsle.vv
 * Integer compare vmsle (vv form)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmsle.vv e8 equal: wrong mask result
 *     2 = vmsle.vv e8 equal: witness changed
 *     3 = vmsle.vv e8 equal: CSR side-effect
 *     4 = vmsle.vv e8 less: wrong mask result
 *     5 = vmsle.vv e8 less: witness changed
 *     6 = vmsle.vv e8 less: CSR side-effect
 *     7 = vmsle.vv e8 greater: wrong mask result
 *     8 = vmsle.vv e8 greater: witness changed
 *     9 = vmsle.vv e8 greater: CSR side-effect
 *    10 = vmsle.vv e8 mixed: wrong mask result
 *    11 = vmsle.vv e8 mixed: witness changed
 *    12 = vmsle.vv e8 mixed: CSR side-effect
 *    13 = vmsle.vv e16 equal: wrong mask result
 *    14 = vmsle.vv e16 equal: witness changed
 *    15 = vmsle.vv e16 equal: CSR side-effect
 *    16 = vmsle.vv e16 less: wrong mask result
 *    17 = vmsle.vv e16 less: witness changed
 *    18 = vmsle.vv e16 less: CSR side-effect
 *    19 = vmsle.vv e16 greater: wrong mask result
 *    20 = vmsle.vv e16 greater: witness changed
 *    21 = vmsle.vv e16 greater: CSR side-effect
 *    22 = vmsle.vv e16 mixed: wrong mask result
 *    23 = vmsle.vv e16 mixed: witness changed
 *    24 = vmsle.vv e16 mixed: CSR side-effect
 *    25 = vmsle.vv e32 equal: wrong mask result
 *    26 = vmsle.vv e32 equal: witness changed
 *    27 = vmsle.vv e32 equal: CSR side-effect
 *    28 = vmsle.vv e32 less: wrong mask result
 *    29 = vmsle.vv e32 less: witness changed
 *    30 = vmsle.vv e32 less: CSR side-effect
 *    31 = vmsle.vv e32 greater: wrong mask result
 *    32 = vmsle.vv e32 greater: witness changed
 *    33 = vmsle.vv e32 greater: CSR side-effect
 *    34 = vmsle.vv e32 mixed: wrong mask result
 *    35 = vmsle.vv e32 mixed: witness changed
 *    36 = vmsle.vv e32 mixed: CSR side-effect
 *    37 = vmsle.vv e32 mixed masked: wrong mask result
 *    38 = vmsle.vv e32 mixed masked: witness changed
 *    39 = vmsle.vv e32 mixed masked: CSR side-effect
 *    40 = vmsle.vv e64 equal: wrong mask result
 *    41 = vmsle.vv e64 equal: witness changed
 *    42 = vmsle.vv e64 equal: CSR side-effect
 *    43 = vmsle.vv e64 less: wrong mask result
 *    44 = vmsle.vv e64 less: witness changed
 *    45 = vmsle.vv e64 less: CSR side-effect
 *    46 = vmsle.vv e64 greater: wrong mask result
 *    47 = vmsle.vv e64 greater: witness changed
 *    48 = vmsle.vv e64 greater: CSR side-effect
 *    49 = vmsle.vv e64 mixed: wrong mask result
 *    50 = vmsle.vv e64 mixed: witness changed
 *    51 = vmsle.vv e64 mixed: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vmsle.vv SEW=8 equal */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    la t1, tc1_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 2
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 4
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED

    /* Test 4-6: vmsle.vv SEW=8 less */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    la t1, tc4_s1
    vle8.v v20, (t1)
    la t1, tc4_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 5
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 4
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-9: vmsle.vv SEW=8 greater */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_s2
    vle8.v v16, (t1)
    la t1, tc7_s1
    vle8.v v20, (t1)
    la t1, tc7_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 8
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 4
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED

    /* Test 10-12: vmsle.vv SEW=8 mixed */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc10_s2
    vle8.v v16, (t1)
    la t1, tc10_s1
    vle8.v v20, (t1)
    la t1, tc10_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 10
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 6
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 11
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 4
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13-15: vmsle.vv SEW=16 equal */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc13_s2
    vle16.v v16, (t1)
    la t1, tc13_s1
    vle16.v v20, (t1)
    la t1, tc13_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 13
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 14
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 8
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED

    /* Test 16-18: vmsle.vv SEW=16 less */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc16_s2
    vle16.v v16, (t1)
    la t1, tc16_s1
    vle16.v v20, (t1)
    la t1, tc16_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 16
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 17
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc16_w, 8
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    /* Test 19-21: vmsle.vv SEW=16 greater */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc19_s2
    vle16.v v16, (t1)
    la t1, tc19_s1
    vle16.v v20, (t1)
    la t1, tc19_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 19
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 20
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc19_w, 8
    SET_TEST_NUM 21
    CHECK_CSRS_UNCHANGED

    /* Test 22-24: vmsle.vv SEW=16 mixed */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc22_s2
    vle16.v v16, (t1)
    la t1, tc22_s1
    vle16.v v20, (t1)
    la t1, tc22_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 22
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 6
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 23
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc22_w, 8
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    /* Test 25-27: vmsle.vv SEW=32 equal */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc25_s2
    vle32.v v16, (t1)
    la t1, tc25_s1
    vle32.v v20, (t1)
    la t1, tc25_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 25
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 26
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc25_w, 16
    SET_TEST_NUM 27
    CHECK_CSRS_UNCHANGED

    /* Test 28-30: vmsle.vv SEW=32 less */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc28_s2
    vle32.v v16, (t1)
    la t1, tc28_s1
    vle32.v v20, (t1)
    la t1, tc28_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 28
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 29
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc28_w, 16
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED

    /* Test 31-33: vmsle.vv SEW=32 greater */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc31_s2
    vle32.v v16, (t1)
    la t1, tc31_s1
    vle32.v v20, (t1)
    la t1, tc31_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 31
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 32
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc31_w, 16
    SET_TEST_NUM 33
    CHECK_CSRS_UNCHANGED

    /* Test 34-36: vmsle.vv SEW=32 mixed */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc34_s2
    vle32.v v16, (t1)
    la t1, tc34_s1
    vle32.v v20, (t1)
    la t1, tc34_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 34
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 6
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 35
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc34_w, 16
    SET_TEST_NUM 36
    CHECK_CSRS_UNCHANGED

    /* Test 37-39: vmsle.vv SEW=32 mixed (masked) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc37_s2
    vle32.v v16, (t1)
    la t1, tc37_s1
    vle32.v v20, (t1)
    la t1, tc37_vdinit
    vlm.v v8, (t1)
    la t1, tc37_mask
    vlm.v v0, (t1)
    la t1, tc37_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20, v0.t
    SET_TEST_NUM 37
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 13
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 38
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc37_w, 16
    SET_TEST_NUM 39
    CHECK_CSRS_UNCHANGED

    /* Test 40-42: vmsle.vv SEW=64 equal */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc40_s2
    vle64.v v16, (t1)
    la t1, tc40_s1
    vle64.v v20, (t1)
    la t1, tc40_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 40
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 41
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc40_w, 32
    SET_TEST_NUM 42
    CHECK_CSRS_UNCHANGED

    /* Test 43-45: vmsle.vv SEW=64 less */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc43_s2
    vle64.v v16, (t1)
    la t1, tc43_s1
    vle64.v v20, (t1)
    la t1, tc43_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 43
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 44
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc43_w, 32
    SET_TEST_NUM 45
    CHECK_CSRS_UNCHANGED

    /* Test 46-48: vmsle.vv SEW=64 greater */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc46_s2
    vle64.v v16, (t1)
    la t1, tc46_s1
    vle64.v v20, (t1)
    la t1, tc46_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 46
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 47
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc46_w, 32
    SET_TEST_NUM 48
    CHECK_CSRS_UNCHANGED

    /* Test 49-51: vmsle.vv SEW=64 mixed */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc49_s2
    vle64.v v16, (t1)
    la t1, tc49_s1
    vle64.v v20, (t1)
    la t1, tc49_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vmsle.vv v8, v16, v20
    SET_TEST_NUM 49
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 6
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 50
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc49_w, 32
    SET_TEST_NUM 51
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_s1:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc4_s2:
    .byte 0x00, 0x01, 0x02, 0x03
tc4_s1:
    .byte 0x01, 0x02, 0x03, 0x04
tc4_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc7_s2:
    .byte 0x02, 0x03, 0x04, 0x05
tc7_s1:
    .byte 0x01, 0x02, 0x03, 0x04
tc7_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc10_s2:
    .byte 0x00, 0xff, 0x80, 0x7f
tc10_s1:
    .byte 0xff, 0x00, 0x7f, 0x80
tc10_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc13_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc13_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc13_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc16_s2:
    .half 0x0000, 0x0001, 0x0002, 0x0003
tc16_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc16_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc19_s2:
    .half 0x0002, 0x0003, 0x0004, 0x0005
tc19_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc19_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc22_s2:
    .half 0x0000, 0xffff, 0x8000, 0x7fff
tc22_s1:
    .half 0xffff, 0x0000, 0x7fff, 0x8000
tc22_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 2
tc25_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc25_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc25_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc28_s2:
    .word 0x00000000, 0x00000001, 0x00000002, 0x00000003
tc28_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc28_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc31_s2:
    .word 0x00000002, 0x00000003, 0x00000004, 0x00000005
tc31_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc31_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc34_s2:
    .word 0x00000000, 0xffffffff, 0x80000000, 0x7fffffff
tc34_s1:
    .word 0xffffffff, 0x00000000, 0x7fffffff, 0x80000000
tc34_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 1
tc37_mask:
    .byte 10
tc37_vdinit:
    .byte 5
.align 2
tc37_s2:
    .word 0x00000005, 0x00000005, 0x00000001, 0x00000001
tc37_s1:
    .word 0x00000005, 0x00000001, 0x00000005, 0x00000005
tc37_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc40_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc40_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc40_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc43_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0x0000000000000002, 0x0000000000000003
tc43_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc43_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc46_s2:
    .dword 0x0000000000000002, 0x0000000000000003, 0x0000000000000004, 0x0000000000000005
tc46_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc46_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc49_s2:
    .dword 0x0000000000000000, 0xffffffffffffffff, 0x8000000000000000, 0x7fffffffffffffff
tc49_s1:
    .dword 0xffffffffffffffff, 0x0000000000000000, 0x7fffffffffffffff, 0x8000000000000000
tc49_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

