/* Auto-generated test for vse32.v
 * Unit-stride store of 32-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vse32.v: store and scalar-verify
 *     2 = vse32.v masked: only active elements stored
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_data
    vle32.v v8, (t1)
    SAVE_CSRS
    la t1, result_buf
    vse32.v v8, (t1)
    SET_TEST_NUM 1
    lwu t2, 0(t1)
    li t3, 1
    FAIL_IF_NE t2, t3
    lwu t2, 4(t1)
    li t3, 2
    FAIL_IF_NE t2, t3
    lwu t2, 8(t1)
    li t3, 3
    FAIL_IF_NE t2, t3
    lwu t2, 12(t1)
    li t3, 4
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    /* Masked store: mask=0b0101, only active written */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc2_src
    vle32.v v8, (t1)
    la t1, tc2_mask
    vlm.v v0, (t1)
    la t1, tc2_buf_init
    vle32.v v16, (t1)
    la t1, result_buf
    vse32.v v16, (t1)
    SAVE_CSRS
    la t1, result_buf
    vse32.v v8, (t1), v0.t
    SET_TEST_NUM 2
    CHECK_MEM result_buf, tc2_exp, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_data:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 1
tc2_mask:
    .byte 5
.align 2
tc2_src:
    .word 0x00000011, 0x00000022, 0x00000033, 0x00000044
tc2_buf_init:
    .word 0x000000ff, 0x000000ff, 0x000000ff, 0x000000ff
tc2_exp:
    .word 0x00000011, 0x000000ff, 0x00000033, 0x000000ff

.align 4
result_buf:  .space 256
witness_buf: .space 256

