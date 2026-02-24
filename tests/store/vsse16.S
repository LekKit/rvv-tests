/* Auto-generated test for vsse16.v
 * Strided store of 16-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsse16.v stride=4: result
 *     2 = vsse16.v stride=4 masked: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_data
    vle16.v v8, (t1)
    la t1, result_buf
    li t2, 4
    SAVE_CSRS
    vsse16.v v8, (t1), t2
    SET_TEST_NUM 1
    vlse16.v v16, (t1), t2
    la t1, witness_buf
    vse16.v v16, (t1)
    CHECK_MEM witness_buf, tc1_data, 8
    CHECK_CSRS_UNCHANGED

    /* Masked strided store: mask=0b0101 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_src
    vle16.v v8, (t1)
    la t1, tc2_mask
    vlm.v v0, (t1)
    vmv.v.i v16, 0
    la t1, result_buf
    li t2, 4
    vsse16.v v16, (t1), t2
    SAVE_CSRS
    la t1, result_buf
    vsse16.v v8, (t1), t2, v0.t
    SET_TEST_NUM 2
    vlse16.v v16, (t1), t2
    la t1, witness_buf
    vse16.v v16, (t1)
    CHECK_MEM witness_buf, tc2_exp, 8
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc2_mask:
    .byte 5
.align 1
tc2_src:
    .half 0x0011, 0x0022, 0x0033, 0x0044
tc2_exp:
    .half 0x0011, 0x0000, 0x0033, 0x0000

.align 4
result_buf:  .space 256
witness_buf: .space 256

