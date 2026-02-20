/* Auto-generated test for vsoxei16.v
 * Indexed store vsoxei16.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsoxei16.v identity: result
 *     2 = vsoxei16.v reverse: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_data
    vle16.v v8, (t1)
    la t1, tc1_idx
    vle16.v v20, (t1)
    la t1, result_buf
    SAVE_CSRS
    vsoxei16.v v8, (t1), v20
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_data
    vle16.v v8, (t1)
    la t1, tc2_idx
    vle16.v v20, (t1)
    la t1, result_buf
    SAVE_CSRS
    vsoxei16.v v8, (t1), v20
    SET_TEST_NUM 2
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc1_idx:
    .half 0x0000, 0x0002, 0x0004, 0x0006
tc1_exp:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc2_data:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc2_idx:
    .half 0x0006, 0x0004, 0x0002, 0x0000
tc2_exp:
    .half 0x0028, 0x001e, 0x0014, 0x000a

.align 4
result_buf:  .space 256
witness_buf: .space 256

