/* Auto-generated test for vle16.v
 * Unit-stride load of 16-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vle16.v: basic load-store round-trip
 *     2 = vle16.v: scalar verify element 0
 *     3 = vle16.v masked: active+inactive elements
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_data
    vle16.v v8, (t1)
    SAVE_CSRS
    SET_TEST_NUM 1
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 8
    CHECK_CSRS_UNCHANGED

    SET_TEST_NUM 2
    la t1, result_buf
    lhu t2, 0(t1)
    li t3, 1
    FAIL_IF_NE t2, t3

    /* Masked load: mask=0b0101, tu/mu policy */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc3_vd
    vle16.v v8, (t1)
    la t1, tc3_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    la t1, tc3_mem
    vle16.v v8, (t1), v0.t
    SET_TEST_NUM 3
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .half 0x0001, 0x0002, 0x0003, 0x0004
.align 1
tc3_mask:
    .byte 5
.align 1
tc3_vd:
    .half 0x00dd, 0x00dd, 0x00dd, 0x00dd
tc3_mem:
    .half 0x00aa, 0x00bb, 0x00cc, 0x00dd
tc3_exp:
    .half 0x00aa, 0x00dd, 0x00cc, 0x00dd

.align 4
result_buf:  .space 256
witness_buf: .space 256

