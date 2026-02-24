/* Auto-generated test for vse16.v
 * Unit-stride store of 16-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vse16.v: store and scalar-verify
 *     2 = vse16.v masked: only active elements stored
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_data
    vle16.v v8, (t1)
    SAVE_CSRS
    la t1, result_buf
    vse16.v v8, (t1)
    SET_TEST_NUM 1
    lhu t2, 0(t1)
    li t3, 1
    FAIL_IF_NE t2, t3
    lhu t2, 2(t1)
    li t3, 2
    FAIL_IF_NE t2, t3
    lhu t2, 4(t1)
    li t3, 3
    FAIL_IF_NE t2, t3
    lhu t2, 6(t1)
    li t3, 4
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    /* Masked store: mask=0b0101, only active written */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_src
    vle16.v v8, (t1)
    la t1, tc2_mask
    vlm.v v0, (t1)
    la t1, tc2_buf_init
    vle16.v v16, (t1)
    la t1, result_buf
    vse16.v v16, (t1)
    SAVE_CSRS
    la t1, result_buf
    vse16.v v8, (t1), v0.t
    SET_TEST_NUM 2
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc2_mask:
    .byte 5
.align 1
tc2_src:
    .half 0x0011, 0x0022, 0x0033, 0x0044
tc2_buf_init:
    .half 0x00ff, 0x00ff, 0x00ff, 0x00ff
tc2_exp:
    .half 0x0011, 0x00ff, 0x0033, 0x00ff

.align 4
result_buf:  .space 256
witness_buf: .space 256

