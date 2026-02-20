/* Auto-generated test for vmsgt.vi
 * Integer compare vmsgt (vi form)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmsgt.vi e8 imm0: wrong mask
 *     2 = vmsgt.vi e8 imm0: CSR
 *     3 = vmsgt.vi e8 imm1: wrong mask
 *     4 = vmsgt.vi e8 imm1: CSR
 *     5 = vmsgt.vi e8 imm_n1: wrong mask
 *     6 = vmsgt.vi e8 imm_n1: CSR
 *     7 = vmsgt.vi e8 imm15: wrong mask
 *     8 = vmsgt.vi e8 imm15: CSR
 *     9 = vmsgt.vi e8 imm_n16: wrong mask
 *    10 = vmsgt.vi e8 imm_n16: CSR
 *    11 = vmsgt.vi e16 imm0: wrong mask
 *    12 = vmsgt.vi e16 imm0: CSR
 *    13 = vmsgt.vi e16 imm1: wrong mask
 *    14 = vmsgt.vi e16 imm1: CSR
 *    15 = vmsgt.vi e16 imm_n1: wrong mask
 *    16 = vmsgt.vi e16 imm_n1: CSR
 *    17 = vmsgt.vi e16 imm15: wrong mask
 *    18 = vmsgt.vi e16 imm15: CSR
 *    19 = vmsgt.vi e16 imm_n16: wrong mask
 *    20 = vmsgt.vi e16 imm_n16: CSR
 *    21 = vmsgt.vi e32 imm0: wrong mask
 *    22 = vmsgt.vi e32 imm0: CSR
 *    23 = vmsgt.vi e32 imm1: wrong mask
 *    24 = vmsgt.vi e32 imm1: CSR
 *    25 = vmsgt.vi e32 imm_n1: wrong mask
 *    26 = vmsgt.vi e32 imm_n1: CSR
 *    27 = vmsgt.vi e32 imm15: wrong mask
 *    28 = vmsgt.vi e32 imm15: CSR
 *    29 = vmsgt.vi e32 imm_n16: wrong mask
 *    30 = vmsgt.vi e32 imm_n16: CSR
 *    31 = vmsgt.vi e64 imm0: wrong mask
 *    32 = vmsgt.vi e64 imm0: CSR
 *    33 = vmsgt.vi e64 imm1: wrong mask
 *    34 = vmsgt.vi e64 imm1: CSR
 *    35 = vmsgt.vi e64 imm_n1: wrong mask
 *    36 = vmsgt.vi e64 imm_n1: CSR
 *    37 = vmsgt.vi e64 imm15: wrong mask
 *    38 = vmsgt.vi e64 imm15: CSR
 *    39 = vmsgt.vi e64 imm_n16: wrong mask
 *    40 = vmsgt.vi e64 imm_n16: CSR
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vmsgt.vi SEW=8 imm0 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 0
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    /* Test 3: vmsgt.vi SEW=8 imm1 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 1
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 14
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED

    /* Test 5: vmsgt.vi SEW=8 imm_n1 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, -1
    SET_TEST_NUM 5
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7: vmsgt.vi SEW=8 imm15 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 15
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED

    /* Test 9: vmsgt.vi SEW=8 imm_n16 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc9_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, -16
    SET_TEST_NUM 9
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED

    /* Test 11: vmsgt.vi SEW=16 imm0 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 0
    SET_TEST_NUM 11
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13: vmsgt.vi SEW=16 imm1 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc13_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 1
    SET_TEST_NUM 13
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 14
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 14
    CHECK_CSRS_UNCHANGED

    /* Test 15: vmsgt.vi SEW=16 imm_n1 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc15_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, -1
    SET_TEST_NUM 15
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 16
    CHECK_CSRS_UNCHANGED

    /* Test 17: vmsgt.vi SEW=16 imm15 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc17_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 15
    SET_TEST_NUM 17
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    /* Test 19: vmsgt.vi SEW=16 imm_n16 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc19_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, -16
    SET_TEST_NUM 19
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 20
    CHECK_CSRS_UNCHANGED

    /* Test 21: vmsgt.vi SEW=32 imm0 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc21_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 0
    SET_TEST_NUM 21
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 22
    CHECK_CSRS_UNCHANGED

    /* Test 23: vmsgt.vi SEW=32 imm1 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc23_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 1
    SET_TEST_NUM 23
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 14
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    /* Test 25: vmsgt.vi SEW=32 imm_n1 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc25_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, -1
    SET_TEST_NUM 25
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 26
    CHECK_CSRS_UNCHANGED

    /* Test 27: vmsgt.vi SEW=32 imm15 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc27_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 15
    SET_TEST_NUM 27
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 28
    CHECK_CSRS_UNCHANGED

    /* Test 29: vmsgt.vi SEW=32 imm_n16 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc29_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, -16
    SET_TEST_NUM 29
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED

    /* Test 31: vmsgt.vi SEW=64 imm0 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc31_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 0
    SET_TEST_NUM 31
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 32
    CHECK_CSRS_UNCHANGED

    /* Test 33: vmsgt.vi SEW=64 imm1 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc33_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 1
    SET_TEST_NUM 33
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 14
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 34
    CHECK_CSRS_UNCHANGED

    /* Test 35: vmsgt.vi SEW=64 imm_n1 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc35_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, -1
    SET_TEST_NUM 35
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 36
    CHECK_CSRS_UNCHANGED

    /* Test 37: vmsgt.vi SEW=64 imm15 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc37_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, 15
    SET_TEST_NUM 37
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 38
    CHECK_CSRS_UNCHANGED

    /* Test 39: vmsgt.vi SEW=64 imm_n16 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc39_s2
    vle64.v v16, (t1)
    SAVE_CSRS
    vmsgt.vi v8, v16, -16
    SET_TEST_NUM 39
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
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
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc9_s2:
    .byte 0x00, 0xff, 0x01, 0xfe
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
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc19_s2:
    .half 0x0000, 0xffff, 0x0001, 0xfffe
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
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc29_s2:
    .word 0x00000000, 0xffffffff, 0x00000001, 0xfffffffe
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
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc39_s2:
    .dword 0x0000000000000000, 0xffffffffffffffff, 0x0000000000000001, 0xfffffffffffffffe

.align 4
result_buf:  .space 256
witness_buf: .space 256

