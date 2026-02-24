/* Auto-generated test for vse64.v
 * Unit-stride store of 64-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vse64.v: store and scalar-verify
 *     2 = vse64.v masked: only active elements stored
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_data
    vle64.v v8, (t1)
    SAVE_CSRS
    la t1, result_buf
    vse64.v v8, (t1)
    SET_TEST_NUM 1
    ld t2, 0(t1)
    li t3, 1
    FAIL_IF_NE t2, t3
    ld t2, 8(t1)
    li t3, 2
    FAIL_IF_NE t2, t3
    ld t2, 16(t1)
    li t3, 3
    FAIL_IF_NE t2, t3
    ld t2, 24(t1)
    li t3, 4
    FAIL_IF_NE t2, t3
    CHECK_CSRS_UNCHANGED

    /* Masked store: mask=0b0101, only active written */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc2_src
    vle64.v v8, (t1)
    la t1, tc2_mask
    vlm.v v0, (t1)
    la t1, tc2_buf_init
    vle64.v v16, (t1)
    la t1, result_buf
    vse64.v v16, (t1)
    SAVE_CSRS
    la t1, result_buf
    vse64.v v8, (t1), v0.t
    SET_TEST_NUM 2
    CHECK_MEM result_buf, tc2_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_data:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 1
tc2_mask:
    .byte 5
.align 3
tc2_src:
    .dword 0x0000000000000011, 0x0000000000000022, 0x0000000000000033, 0x0000000000000044
tc2_buf_init:
    .dword 0x00000000000000ff, 0x00000000000000ff, 0x00000000000000ff, 0x00000000000000ff
tc2_exp:
    .dword 0x0000000000000011, 0x00000000000000ff, 0x0000000000000033, 0x00000000000000ff

.align 4
result_buf:  .space 256
witness_buf: .space 256

