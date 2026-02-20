/* Auto-generated test for vluxei8.v
 * Indexed load vluxei8.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vluxei8.v identity: result
 *     2 = vluxei8.v reverse: result
 *     3 = vluxei8.v splat: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_idx
    vle8.v v20, (t1)
    la t1, tc1_data
    SAVE_CSRS
    vluxei8.v v8, (t1), v20
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
    vluxei8.v v8, (t1), v20
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
    vluxei8.v v8, (t1), v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
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

.align 4
result_buf:  .space 256
witness_buf: .space 256

