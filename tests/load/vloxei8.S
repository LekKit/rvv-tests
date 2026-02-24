/* Auto-generated test for vloxei8.v
 * Indexed load vloxei8.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vloxei8.v identity: result
 *     2 = vloxei8.v reverse: result
 *     3 = vloxei8.v splat: result
 *     4 = vloxei8.v masked: active+inactive
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_idx
    vle8.v v20, (t1)
    la t1, tc1_data
    SAVE_CSRS
    vloxei8.v v8, (t1), v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_idx
    vle8.v v20, (t1)
    la t1, tc2_data
    SAVE_CSRS
    vloxei8.v v8, (t1), v20
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_idx
    vle8.v v20, (t1)
    la t1, tc3_data
    SAVE_CSRS
    vloxei8.v v8, (t1), v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED

    /* Masked indexed load: mask=0b0101 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_vd
    vle8.v v8, (t1)
    la t1, tc4_idx
    vle8.v v20, (t1)
    la t1, tc4_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    la t1, tc4_data
    vloxei8.v v8, (t1), v20, v0.t
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .byte 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11
tc1_idx:
    .byte 0x00, 0x01, 0x02, 0x03
tc1_exp:
    .byte 0x0a, 0x0b, 0x0c, 0x0d
.align 1
tc2_data:
    .byte 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11
tc2_idx:
    .byte 0x03, 0x02, 0x01, 0x00
tc2_exp:
    .byte 0x0d, 0x0c, 0x0b, 0x0a
.align 1
tc3_data:
    .byte 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11
tc3_idx:
    .byte 0x02, 0x02, 0x02, 0x02
tc3_exp:
    .byte 0x0c, 0x0c, 0x0c, 0x0c
.align 1
tc4_mask:
    .byte 5
.align 1
tc4_vd:
    .byte 0xdd, 0xdd, 0xdd, 0xdd
tc4_data:
    .byte 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11
tc4_idx:
    .byte 0x00, 0x01, 0x02, 0x03
tc4_exp:
    .byte 0x0a, 0xdd, 0x0c, 0xdd

.align 4
result_buf:  .space 256
witness_buf: .space 256

