/* Auto-generated test for vsoxei8.v
 * Indexed store vsoxei8.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsoxei8.v identity: result
 *     2 = vsoxei8.v reverse: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_data
    vle8.v v8, (t1)
    la t1, tc1_idx
    vle8.v v20, (t1)
    la t1, result_buf
    SAVE_CSRS
    vsoxei8.v v8, (t1), v20
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_data
    vle8.v v8, (t1)
    la t1, tc2_idx
    vle8.v v20, (t1)
    la t1, result_buf
    SAVE_CSRS
    vsoxei8.v v8, (t1), v20
    SET_TEST_NUM 2
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_idx:
    .byte 0x00, 0x01, 0x02, 0x03
tc1_exp:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc2_data:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc2_idx:
    .byte 0x03, 0x02, 0x01, 0x00
tc2_exp:
    .byte 0x28, 0x1e, 0x14, 0x0a

.align 4
result_buf:  .space 256
witness_buf: .space 256

