/* Auto-generated test for vloxei32.v
 * Indexed load vloxei32.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vloxei32.v identity: result
 *     2 = vloxei32.v reverse: result
 *     3 = vloxei32.v splat: result
 *     4 = vloxei32.v masked: active+inactive
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_idx
    vle32.v v20, (t1)
    la t1, tc1_data
    SAVE_CSRS
    vloxei32.v v8, (t1), v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc2_idx
    vle32.v v20, (t1)
    la t1, tc2_data
    SAVE_CSRS
    vloxei32.v v8, (t1), v20
    SET_TEST_NUM 2
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_idx
    vle32.v v20, (t1)
    la t1, tc3_data
    SAVE_CSRS
    vloxei32.v v8, (t1), v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 16
    CHECK_CSRS_UNCHANGED

    /* Masked indexed load: mask=0b0101 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc4_vd
    vle32.v v8, (t1)
    la t1, tc4_idx
    vle32.v v20, (t1)
    la t1, tc4_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    la t1, tc4_data
    vloxei32.v v8, (t1), v20, v0.t
    SET_TEST_NUM 4
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_data:
    .word 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010, 0x00000011
tc1_idx:
    .word 0x00000000, 0x00000004, 0x00000008, 0x0000000c
tc1_exp:
    .word 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d
.align 2
tc2_data:
    .word 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010, 0x00000011
tc2_idx:
    .word 0x0000000c, 0x00000008, 0x00000004, 0x00000000
tc2_exp:
    .word 0x0000000d, 0x0000000c, 0x0000000b, 0x0000000a
.align 2
tc3_data:
    .word 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010, 0x00000011
tc3_idx:
    .word 0x00000008, 0x00000008, 0x00000008, 0x00000008
tc3_exp:
    .word 0x0000000c, 0x0000000c, 0x0000000c, 0x0000000c
.align 1
tc4_mask:
    .byte 5
.align 2
tc4_vd:
    .word 0x000000dd, 0x000000dd, 0x000000dd, 0x000000dd
tc4_data:
    .word 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010, 0x00000011
tc4_idx:
    .word 0x00000000, 0x00000004, 0x00000008, 0x0000000c
tc4_exp:
    .word 0x0000000a, 0x000000dd, 0x0000000c, 0x000000dd

.align 4
result_buf:  .space 256
witness_buf: .space 256

