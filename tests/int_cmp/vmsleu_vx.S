/* Auto-generated test for vmsleu.vx
 * Integer compare vmsleu (vx form)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmsleu.vx e8 basic: wrong mask
 *     2 = vmsleu.vx e8 basic: CSR
 *     3 = vmsleu.vx e8 zero_sc: wrong mask
 *     4 = vmsleu.vx e8 zero_sc: CSR
 *     5 = vmsleu.vx e8 max_sc: wrong mask
 *     6 = vmsleu.vx e8 max_sc: CSR
 *     7 = vmsleu.vx e8 max_vec: wrong mask
 *     8 = vmsleu.vx e8 max_vec: CSR
 *     9 = vmsleu.vx e8 signed_sc: wrong mask
 *    10 = vmsleu.vx e8 signed_sc: CSR
 *    11 = vmsleu.vx e16 basic: wrong mask
 *    12 = vmsleu.vx e16 basic: CSR
 *    13 = vmsleu.vx e16 zero_sc: wrong mask
 *    14 = vmsleu.vx e16 zero_sc: CSR
 *    15 = vmsleu.vx e16 max_sc: wrong mask
 *    16 = vmsleu.vx e16 max_sc: CSR
 *    17 = vmsleu.vx e16 max_vec: wrong mask
 *    18 = vmsleu.vx e16 max_vec: CSR
 *    19 = vmsleu.vx e16 signed_sc: wrong mask
 *    20 = vmsleu.vx e16 signed_sc: CSR
 *    21 = vmsleu.vx e32 basic: wrong mask
 *    22 = vmsleu.vx e32 basic: CSR
 *    23 = vmsleu.vx e32 zero_sc: wrong mask
 *    24 = vmsleu.vx e32 zero_sc: CSR
 *    25 = vmsleu.vx e32 max_sc: wrong mask
 *    26 = vmsleu.vx e32 max_sc: CSR
 *    27 = vmsleu.vx e32 max_vec: wrong mask
 *    28 = vmsleu.vx e32 max_vec: CSR
 *    29 = vmsleu.vx e32 signed_sc: wrong mask
 *    30 = vmsleu.vx e32 signed_sc: CSR
 *    31 = vmsleu.vx e64 basic: wrong mask
 *    32 = vmsleu.vx e64 basic: CSR
 *    33 = vmsleu.vx e64 zero_sc: wrong mask
 *    34 = vmsleu.vx e64 zero_sc: CSR
 *    35 = vmsleu.vx e64 max_sc: wrong mask
 *    36 = vmsleu.vx e64 max_sc: CSR
 *    37 = vmsleu.vx e64 max_vec: wrong mask
 *    38 = vmsleu.vx e64 max_vec: CSR
 *    39 = vmsleu.vx e64 signed_sc: wrong mask
 *    40 = vmsleu.vx e64 signed_sc: CSR
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vmsleu.vx SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    li a0, 0x0a
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    /* Test 3: vmsleu.vx SEW=8 zero_sc */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    li a0, 0x00
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED

    /* Test 5: vmsleu.vx SEW=8 max_sc */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_s2
    vle8.v v16, (t1)
    li a0, 0xff
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 5
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7: vmsleu.vx SEW=8 max_vec */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_s2
    vle8.v v16, (t1)
    li a0, 0x01
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED

    /* Test 9: vmsleu.vx SEW=8 signed_sc */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc9_s2
    vle8.v v16, (t1)
    li a0, 0x80
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 9
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 7
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED

    /* Test 11: vmsleu.vx SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_s2
    vle16.v v16, (t1)
    li a0, 0x000a
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 11
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13: vmsleu.vx SEW=16 zero_sc */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc13_s2
    vle16.v v16, (t1)
    li a0, 0x0000
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 13
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 14
    CHECK_CSRS_UNCHANGED

    /* Test 15: vmsleu.vx SEW=16 max_sc */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc15_s2
    vle16.v v16, (t1)
    li a0, 0xffff
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 15
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 16
    CHECK_CSRS_UNCHANGED

    /* Test 17: vmsleu.vx SEW=16 max_vec */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc17_s2
    vle16.v v16, (t1)
    li a0, 0x0001
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 17
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    /* Test 19: vmsleu.vx SEW=16 signed_sc */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc19_s2
    vle16.v v16, (t1)
    li a0, 0x8000
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 19
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 7
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 20
    CHECK_CSRS_UNCHANGED

    /* Test 21: vmsleu.vx SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc21_s2
    vle32.v v16, (t1)
    li a0, 0x0000000a
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 21
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 22
    CHECK_CSRS_UNCHANGED

    /* Test 23: vmsleu.vx SEW=32 zero_sc */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc23_s2
    vle32.v v16, (t1)
    li a0, 0x00000000
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 23
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    /* Test 25: vmsleu.vx SEW=32 max_sc */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc25_s2
    vle32.v v16, (t1)
    li a0, 0xffffffff
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 25
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 26
    CHECK_CSRS_UNCHANGED

    /* Test 27: vmsleu.vx SEW=32 max_vec */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc27_s2
    vle32.v v16, (t1)
    li a0, 0x00000001
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 27
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 28
    CHECK_CSRS_UNCHANGED

    /* Test 29: vmsleu.vx SEW=32 signed_sc */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc29_s2
    vle32.v v16, (t1)
    li a0, 0x80000000
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 29
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 7
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED

    /* Test 31: vmsleu.vx SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc31_s2
    vle64.v v16, (t1)
    li a0, 0x000000000000000a
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 31
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 32
    CHECK_CSRS_UNCHANGED

    /* Test 33: vmsleu.vx SEW=64 zero_sc */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc33_s2
    vle64.v v16, (t1)
    li a0, 0x0000000000000000
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 33
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 34
    CHECK_CSRS_UNCHANGED

    /* Test 35: vmsleu.vx SEW=64 max_sc */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc35_s2
    vle64.v v16, (t1)
    li a0, 0xffffffffffffffff
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 35
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 36
    CHECK_CSRS_UNCHANGED

    /* Test 37: vmsleu.vx SEW=64 max_vec */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc37_s2
    vle64.v v16, (t1)
    li a0, 0x0000000000000001
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 37
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 38
    CHECK_CSRS_UNCHANGED

    /* Test 39: vmsleu.vx SEW=64 signed_sc */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc39_s2
    vle64.v v16, (t1)
    li a0, 0x8000000000000000
    SAVE_CSRS
    vmsleu.vx v8, v16, a0
    SET_TEST_NUM 39
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 7
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 40
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc3_s2:
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc5_s2:
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc7_s2:
    .byte 0xff, 0xff, 0xff, 0xff
.align 1
tc9_s2:
    .byte 0x00, 0x01, 0x80, 0xff
.align 1
tc11_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc13_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc15_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc17_s2:
    .half 0xffff, 0xffff, 0xffff, 0xffff
.align 1
tc19_s2:
    .half 0x0000, 0x0001, 0x8000, 0xffff
.align 2
tc21_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc23_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc25_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc27_s2:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
.align 2
tc29_s2:
    .word 0x00000000, 0x00000001, 0x80000000, 0xffffffff
.align 3
tc31_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc33_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc35_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc37_s2:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
.align 3
tc39_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0x8000000000000000, 0xffffffffffffffff

.align 4
result_buf:  .space 256
witness_buf: .space 256

