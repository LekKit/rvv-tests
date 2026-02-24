/* Auto-generated test for vsse64.v
 * Strided store of 64-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsse64.v stride=16: result
 *     2 = vsse64.v stride=16 masked: result
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

    /* Masked strided store: mask=0b0101 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc2_src
    vle64.v v8, (t1)
    la t1, tc2_mask
    vlm.v v0, (t1)
    vmv.v.i v16, 0
    la t1, result_buf
    li t2, 16
    vsse64.v v16, (t1), t2
    SAVE_CSRS
    la t1, result_buf
    vsse64.v v8, (t1), t2, v0.t
    SET_TEST_NUM 2
    vlse64.v v16, (t1), t2
    la t1, witness_buf
    vse64.v v16, (t1)
    CHECK_MEM witness_buf, tc2_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_data:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 1
tc2_mask:
    .byte 5
.align 3
tc2_src:
    .dword 0x0000000000000011, 0x0000000000000022, 0x0000000000000033, 0x0000000000000044
tc2_exp:
    .dword 0x0000000000000011, 0x0000000000000000, 0x0000000000000033, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

