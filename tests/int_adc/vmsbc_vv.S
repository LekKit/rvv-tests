/* Auto-generated test for vmsbc.vv
 * Borrow-out detection (no borrow-in)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmsbc.vv e8 basic: result
 *     2 = vmsbc.vv e8 under: result
 *     3 = vmsbc.vv e8 equal: result
 *     4 = vmsbc.vv e16 basic: result
 *     5 = vmsbc.vv e16 under: result
 *     6 = vmsbc.vv e16 equal: result
 *     7 = vmsbc.vv e32 basic: result
 *     8 = vmsbc.vv e32 under: result
 *     9 = vmsbc.vv e32 equal: result
 *    10 = vmsbc.vv e64 basic: result
 *    11 = vmsbc.vv e64 under: result
 *    12 = vmsbc.vv e64 equal: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_s2
    vle8.v v16, (t1)
    la t1, tc2_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 2
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_s2
    vle8.v v16, (t1)
    la t1, tc3_s1
    vle8.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc4_s2
    vle16.v v16, (t1)
    la t1, tc4_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_s2
    vle16.v v16, (t1)
    la t1, tc5_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 5
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_s2
    vle16.v v16, (t1)
    la t1, tc6_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 6
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc7_s2
    vle32.v v16, (t1)
    la t1, tc7_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc8_s2
    vle32.v v16, (t1)
    la t1, tc8_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 8
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_s2
    vle32.v v16, (t1)
    la t1, tc9_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 9
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc10_s2
    vle64.v v16, (t1)
    la t1, tc10_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 10
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc11_s2
    vle64.v v16, (t1)
    la t1, tc11_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 11
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc12_s2
    vle64.v v16, (t1)
    la t1, tc12_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vmsbc.vv v8, v16, v20
    SET_TEST_NUM 12
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_s1:
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc2_s2:
    .byte 0x00, 0x01, 0x00, 0x02
tc2_s1:
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc3_s2:
    .byte 0x05, 0x05, 0x05, 0x05
tc3_s1:
    .byte 0x05, 0x05, 0x05, 0x05
.align 1
tc4_s2:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc4_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc5_s2:
    .half 0x0000, 0x0001, 0x0000, 0x0002
tc5_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc6_s2:
    .half 0x0005, 0x0005, 0x0005, 0x0005
tc6_s1:
    .half 0x0005, 0x0005, 0x0005, 0x0005
.align 2
tc7_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc7_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc8_s2:
    .word 0x00000000, 0x00000001, 0x00000000, 0x00000002
tc8_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 2
tc9_s2:
    .word 0x00000005, 0x00000005, 0x00000005, 0x00000005
tc9_s1:
    .word 0x00000005, 0x00000005, 0x00000005, 0x00000005
.align 3
tc10_s2:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc10_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc11_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0x0000000000000000, 0x0000000000000002
tc11_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc12_s2:
    .dword 0x0000000000000005, 0x0000000000000005, 0x0000000000000005, 0x0000000000000005
tc12_s1:
    .dword 0x0000000000000005, 0x0000000000000005, 0x0000000000000005, 0x0000000000000005

.align 4
result_buf:  .space 256
witness_buf: .space 256

