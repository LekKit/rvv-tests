/* Auto-generated test for vsuxei64.v
 * Indexed store vsuxei64.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsuxei64.v identity: result
 *     2 = vsuxei64.v reverse: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_data
    vle64.v v8, (t1)
    la t1, tc1_idx
    vle64.v v20, (t1)
    la t1, result_buf
    SAVE_CSRS
    vsuxei64.v v8, (t1), v20
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc2_data
    vle64.v v8, (t1)
    la t1, tc2_idx
    vle64.v v20, (t1)
    la t1, result_buf
    SAVE_CSRS
    vsuxei64.v v8, (t1), v20
    SET_TEST_NUM 2
    CHECK_MEM result_buf, tc2_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_data:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc1_idx:
    .dword 0x0000000000000000, 0x0000000000000008, 0x0000000000000010, 0x0000000000000018
tc1_exp:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 3
tc2_data:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc2_idx:
    .dword 0x0000000000000018, 0x0000000000000010, 0x0000000000000008, 0x0000000000000000
tc2_exp:
    .dword 0x0000000000000028, 0x000000000000001e, 0x0000000000000014, 0x000000000000000a

.align 4
result_buf:  .space 256
witness_buf: .space 256

