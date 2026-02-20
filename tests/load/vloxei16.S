/* Auto-generated test for vloxei16.v
 * Indexed load vloxei16.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vloxei16.v identity: result
 *     2 = vloxei16.v reverse: result
 *     3 = vloxei16.v splat: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_idx
    vle16.v v20, (t1)
    la t1, tc1_data
    SAVE_CSRS
    vloxei16.v v8, (t1), v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_idx
    vle16.v v20, (t1)
    la t1, tc2_data
    SAVE_CSRS
    vloxei16.v v8, (t1), v20
    SET_TEST_NUM 2
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc3_idx
    vle16.v v20, (t1)
    la t1, tc3_data
    SAVE_CSRS
    vloxei16.v v8, (t1), v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .half 0x000a, 0x000b, 0x000c, 0x000d, 0x000e, 0x000f, 0x0010, 0x0011
tc1_idx:
    .half 0x0000, 0x0002, 0x0004, 0x0006
tc1_exp:
    .half 0x000a, 0x000b, 0x000c, 0x000d
.align 1
tc2_data:
    .half 0x000a, 0x000b, 0x000c, 0x000d, 0x000e, 0x000f, 0x0010, 0x0011
tc2_idx:
    .half 0x0006, 0x0004, 0x0002, 0x0000
tc2_exp:
    .half 0x000d, 0x000c, 0x000b, 0x000a
.align 1
tc3_data:
    .half 0x000a, 0x000b, 0x000c, 0x000d, 0x000e, 0x000f, 0x0010, 0x0011
tc3_idx:
    .half 0x0004, 0x0004, 0x0004, 0x0004
tc3_exp:
    .half 0x000c, 0x000c, 0x000c, 0x000c

.align 4
result_buf:  .space 256
witness_buf: .space 256

