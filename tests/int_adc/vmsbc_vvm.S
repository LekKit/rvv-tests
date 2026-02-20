/* Auto-generated test for vmsbc.vvm
 * Borrow-out detection with borrow-in
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmsbc.vvm e8 basic: result
 *     2 = vmsbc.vvm e8 carry: result
 *     3 = vmsbc.vvm e8 overflow: result
 *     4 = vmsbc.vvm e8 mixed: result
 *     5 = vmsbc.vvm e16 basic: result
 *     6 = vmsbc.vvm e16 carry: result
 *     7 = vmsbc.vvm e16 overflow: result
 *     8 = vmsbc.vvm e16 mixed: result
 *     9 = vmsbc.vvm e32 basic: result
 *    10 = vmsbc.vvm e32 carry: result
 *    11 = vmsbc.vvm e32 overflow: result
 *    12 = vmsbc.vvm e32 mixed: result
 *    13 = vmsbc.vvm e64 basic: result
 *    14 = vmsbc.vvm e64 carry: result
 *    15 = vmsbc.vvm e64 overflow: result
 *    16 = vmsbc.vvm e64 mixed: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    la t1, tc1_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_s2
    vle8.v v16, (t1)
    la t1, tc2_s1
    vle8.v v20, (t1)
    la t1, tc2_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
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
    la t1, tc3_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 8
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    la t1, tc4_s1
    vle8.v v20, (t1)
    la t1, tc4_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 4
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 5
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_s2
    vle16.v v16, (t1)
    la t1, tc5_s1
    vle16.v v20, (t1)
    la t1, tc5_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
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
    la t1, tc6_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 6
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_s2
    vle16.v v16, (t1)
    la t1, tc7_s1
    vle16.v v20, (t1)
    la t1, tc7_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 8
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc8_s2
    vle16.v v16, (t1)
    la t1, tc8_s1
    vle16.v v20, (t1)
    la t1, tc8_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 8
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 5
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_s2
    vle32.v v16, (t1)
    la t1, tc9_s1
    vle32.v v20, (t1)
    la t1, tc9_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 9
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc10_s2
    vle32.v v16, (t1)
    la t1, tc10_s1
    vle32.v v20, (t1)
    la t1, tc10_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 10
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_s2
    vle32.v v16, (t1)
    la t1, tc11_s1
    vle32.v v20, (t1)
    la t1, tc11_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 11
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 8
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc12_s2
    vle32.v v16, (t1)
    la t1, tc12_s1
    vle32.v v20, (t1)
    la t1, tc12_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 12
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 5
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc13_s2
    vle64.v v16, (t1)
    la t1, tc13_s1
    vle64.v v20, (t1)
    la t1, tc13_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 13
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc14_s2
    vle64.v v16, (t1)
    la t1, tc14_s1
    vle64.v v20, (t1)
    la t1, tc14_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 14
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 15
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc15_s2
    vle64.v v16, (t1)
    la t1, tc15_s1
    vle64.v v20, (t1)
    la t1, tc15_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 15
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 8
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_s2
    vle64.v v16, (t1)
    la t1, tc16_s1
    vle64.v v20, (t1)
    la t1, tc16_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmsbc.vvm v8, v16, v20, v0
    SET_TEST_NUM 16
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 5
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_mask:
    .byte 0
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_s1:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc2_mask:
    .byte 15
.align 1
tc2_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc2_s1:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc3_mask:
    .byte 10
.align 1
tc3_s2:
    .byte 0xff, 0xff, 0xfe, 0x00
tc3_s1:
    .byte 0x01, 0x02, 0x03, 0xff
.align 1
tc4_mask:
    .byte 5
.align 1
tc4_s2:
    .byte 0x00, 0xff, 0x01, 0x7f
tc4_s1:
    .byte 0xff, 0x00, 0xfe, 0x01
.align 1
tc5_mask:
    .byte 0
.align 1
tc5_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc5_s1:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc6_mask:
    .byte 15
.align 1
tc6_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc6_s1:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc7_mask:
    .byte 10
.align 1
tc7_s2:
    .half 0xffff, 0xffff, 0xfffe, 0x0000
tc7_s1:
    .half 0x0001, 0x0002, 0x0003, 0xffff
.align 1
tc8_mask:
    .byte 5
.align 1
tc8_s2:
    .half 0x0000, 0xffff, 0x0001, 0x7fff
tc8_s1:
    .half 0xffff, 0x0000, 0xfffe, 0x0001
.align 1
tc9_mask:
    .byte 0
.align 2
tc9_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc9_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 1
tc10_mask:
    .byte 15
.align 2
tc10_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc10_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 1
tc11_mask:
    .byte 10
.align 2
tc11_s2:
    .word 0xffffffff, 0xffffffff, 0xfffffffe, 0x00000000
tc11_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0xffffffff
.align 1
tc12_mask:
    .byte 5
.align 2
tc12_s2:
    .word 0x00000000, 0xffffffff, 0x00000001, 0x7fffffff
tc12_s1:
    .word 0xffffffff, 0x00000000, 0xfffffffe, 0x00000001
.align 1
tc13_mask:
    .byte 0
.align 3
tc13_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc13_s1:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 1
tc14_mask:
    .byte 15
.align 3
tc14_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc14_s1:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 1
tc15_mask:
    .byte 10
.align 3
tc15_s2:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xfffffffffffffffe, 0x0000000000000000
tc15_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0xffffffffffffffff
.align 1
tc16_mask:
    .byte 5
.align 3
tc16_s2:
    .dword 0x0000000000000000, 0xffffffffffffffff, 0x0000000000000001, 0x7fffffffffffffff
tc16_s1:
    .dword 0xffffffffffffffff, 0x0000000000000000, 0xfffffffffffffffe, 0x0000000000000001

.align 4
result_buf:  .space 256
witness_buf: .space 256

