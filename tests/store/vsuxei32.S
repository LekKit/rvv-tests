/* Auto-generated test for vsuxei32.v
 * Indexed store vsuxei32.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsuxei32.v identity: result
 *     2 = vsuxei32.v reverse: result
 *     3 = vsuxei32.v masked: only active stored
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_data
    vle32.v v8, (t1)
    la t1, tc1_idx
    vle32.v v20, (t1)
    la t1, result_buf
    SAVE_CSRS
    vsuxei32.v v8, (t1), v20
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc2_data
    vle32.v v8, (t1)
    la t1, tc2_idx
    vle32.v v20, (t1)
    la t1, result_buf
    SAVE_CSRS
    vsuxei32.v v8, (t1), v20
    SET_TEST_NUM 2
    CHECK_MEM result_buf, tc2_exp, 16
    CHECK_CSRS_UNCHANGED

    /* Masked indexed store: mask=0b0101 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_src
    vle32.v v8, (t1)
    la t1, tc3_idx
    vle32.v v20, (t1)
    la t1, tc3_mask
    vlm.v v0, (t1)
    la t1, tc3_buf_init
    vle32.v v16, (t1)
    la t1, result_buf
    vse32.v v16, (t1)
    SAVE_CSRS
    la t1, result_buf
    vsuxei32.v v8, (t1), v20, v0.t
    SET_TEST_NUM 3
    CHECK_MEM result_buf, tc3_exp, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_data:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc1_idx:
    .word 0x00000000, 0x00000004, 0x00000008, 0x0000000c
tc1_exp:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 2
tc2_data:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc2_idx:
    .word 0x0000000c, 0x00000008, 0x00000004, 0x00000000
tc2_exp:
    .word 0x00000028, 0x0000001e, 0x00000014, 0x0000000a
.align 1
tc3_mask:
    .byte 5
.align 2
tc3_src:
    .word 0x00000011, 0x00000022, 0x00000033, 0x00000044
tc3_buf_init:
    .word 0x000000ff, 0x000000ff, 0x000000ff, 0x000000ff
tc3_idx:
    .word 0x00000000, 0x00000004, 0x00000008, 0x0000000c
tc3_exp:
    .word 0x00000011, 0x000000ff, 0x00000033, 0x000000ff

.align 4
result_buf:  .space 256
witness_buf: .space 256

