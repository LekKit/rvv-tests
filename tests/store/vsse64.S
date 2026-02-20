/* Auto-generated test for vsse64.v
 * Strided store of 64-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsse64.v stride=16: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_data
    vle64.v v8, (t1)
    la t1, result_buf
    li t2, 16
    SAVE_CSRS
    vsse64.v v8, (t1), t2
    SET_TEST_NUM 1
    vlse64.v v16, (t1), t2
    la t1, witness_buf
    vse64.v v16, (t1)
    CHECK_MEM witness_buf, tc1_data, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_data:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028

.align 4
result_buf:  .space 256
witness_buf: .space 256

