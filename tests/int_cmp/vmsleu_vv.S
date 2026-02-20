/* Auto-generated test for vmsleu.vv
 * Integer compare vmsleu (vv form)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmsleu.vv e8 equal: wrong mask result
 *     2 = vmsleu.vv e8 equal: CSR side-effect
 *     3 = vmsleu.vv e8 less: wrong mask result
 *     4 = vmsleu.vv e8 less: CSR side-effect
 *     5 = vmsleu.vv e8 greater: wrong mask result
 *     6 = vmsleu.vv e8 greater: CSR side-effect
 *     7 = vmsleu.vv e8 mixed: wrong mask result
 *     8 = vmsleu.vv e8 mixed: CSR side-effect
 *     9 = vmsleu.vv e16 equal: wrong mask result
 *    10 = vmsleu.vv e16 equal: CSR side-effect
 *    11 = vmsleu.vv e16 less: wrong mask result
 *    12 = vmsleu.vv e16 less: CSR side-effect
 *    13 = vmsleu.vv e16 greater: wrong mask result
 *    14 = vmsleu.vv e16 greater: CSR side-effect
 *    15 = vmsleu.vv e16 mixed: wrong mask result
 *    16 = vmsleu.vv e16 mixed: CSR side-effect
 *    17 = vmsleu.vv e32 equal: wrong mask result
 *    18 = vmsleu.vv e32 equal: CSR side-effect
 *    19 = vmsleu.vv e32 less: wrong mask result
 *    20 = vmsleu.vv e32 less: CSR side-effect
 *    21 = vmsleu.vv e32 greater: wrong mask result
 *    22 = vmsleu.vv e32 greater: CSR side-effect
 *    23 = vmsleu.vv e32 mixed: wrong mask result
 *    24 = vmsleu.vv e32 mixed: CSR side-effect
 *    25 = vmsleu.vv e64 equal: wrong mask result
 *    26 = vmsleu.vv e64 equal: CSR side-effect
 *    27 = vmsleu.vv e64 less: wrong mask result
 *    28 = vmsleu.vv e64 less: CSR side-effect
 *    29 = vmsleu.vv e64 greater: wrong mask result
 *    30 = vmsleu.vv e64 greater: CSR side-effect
 *    31 = vmsleu.vv e64 mixed: wrong mask result
 *    32 = vmsleu.vv e64 mixed: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vmsleu.vv SEW=8 equal */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    /* Test 3-4: vmsleu.vv SEW=8 less */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    la t1, tc3_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED

    /* Test 5-6: vmsleu.vv SEW=8 greater */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_s2
    vle8.v v16, (t1)
    la t1, tc5_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 5
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-8: vmsleu.vv SEW=8 mixed */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_s2
    vle8.v v16, (t1)
    la t1, tc7_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 9
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED

    /* Test 9-10: vmsleu.vv SEW=16 equal */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_s2
    vle16.v v16, (t1)
    la t1, tc9_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 9
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED

    /* Test 11-12: vmsleu.vv SEW=16 less */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_s2
    vle16.v v16, (t1)
    la t1, tc11_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 11
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13-14: vmsleu.vv SEW=16 greater */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc13_s2
    vle16.v v16, (t1)
    la t1, tc13_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 13
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 14
    CHECK_CSRS_UNCHANGED

    /* Test 15-16: vmsleu.vv SEW=16 mixed */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc15_s2
    vle16.v v16, (t1)
    la t1, tc15_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 15
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 9
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 16
    CHECK_CSRS_UNCHANGED

    /* Test 17-18: vmsleu.vv SEW=32 equal */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc17_s2
    vle32.v v16, (t1)
    la t1, tc17_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 17
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    /* Test 19-20: vmsleu.vv SEW=32 less */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc19_s2
    vle32.v v16, (t1)
    la t1, tc19_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 19
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 20
    CHECK_CSRS_UNCHANGED

    /* Test 21-22: vmsleu.vv SEW=32 greater */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc21_s2
    vle32.v v16, (t1)
    la t1, tc21_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 21
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 22
    CHECK_CSRS_UNCHANGED

    /* Test 23-24: vmsleu.vv SEW=32 mixed */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc23_s2
    vle32.v v16, (t1)
    la t1, tc23_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 23
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 9
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    /* Test 25-26: vmsleu.vv SEW=64 equal */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc25_s2
    vle64.v v16, (t1)
    la t1, tc25_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 25
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 26
    CHECK_CSRS_UNCHANGED

    /* Test 27-28: vmsleu.vv SEW=64 less */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc27_s2
    vle64.v v16, (t1)
    la t1, tc27_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 27
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 28
    CHECK_CSRS_UNCHANGED

    /* Test 29-30: vmsleu.vv SEW=64 greater */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc29_s2
    vle64.v v16, (t1)
    la t1, tc29_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 29
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED

    /* Test 31-32: vmsleu.vv SEW=64 mixed */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc31_s2
    vle64.v v16, (t1)
    la t1, tc31_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vmsleu.vv v8, v16, v20
    SET_TEST_NUM 31
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 9
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_s1:
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc3_s2:
    .byte 0x00, 0x01, 0x02, 0x03
tc3_s1:
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc5_s2:
    .byte 0x02, 0x03, 0x04, 0x05
tc5_s1:
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc7_s2:
    .byte 0x00, 0xff, 0x80, 0x7f
tc7_s1:
    .byte 0xff, 0x00, 0x7f, 0x80
.align 1
tc9_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc9_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc11_s2:
    .half 0x0000, 0x0001, 0x0002, 0x0003
tc11_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc13_s2:
    .half 0x0002, 0x0003, 0x0004, 0x0005
tc13_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc15_s2:
    .half 0x0000, 0xffff, 0x8000, 0x7fff
tc15_s1:
    .half 0xffff, 0x0000, 0x7fff, 0x8000
.align 2
tc17_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc17_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc19_s2:
    .word 0x00000000, 0x00000001, 0x00000002, 0x00000003
tc19_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc21_s2:
    .word 0x00000002, 0x00000003, 0x00000004, 0x00000005
tc21_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc23_s2:
    .word 0x00000000, 0xffffffff, 0x80000000, 0x7fffffff
tc23_s1:
    .word 0xffffffff, 0x00000000, 0x7fffffff, 0x80000000
.align 3
tc25_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc25_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc27_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0x0000000000000002, 0x0000000000000003
tc27_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc29_s2:
    .dword 0x0000000000000002, 0x0000000000000003, 0x0000000000000004, 0x0000000000000005
tc29_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc31_s2:
    .dword 0x0000000000000000, 0xffffffffffffffff, 0x8000000000000000, 0x7fffffffffffffff
tc31_s1:
    .dword 0xffffffffffffffff, 0x0000000000000000, 0x7fffffffffffffff, 0x8000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

