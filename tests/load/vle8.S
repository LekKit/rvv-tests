/* Auto-generated test for vle8.v
 * Unit-stride load of 8-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vle8.v: basic load-store round-trip
 *     2 = vle8.v: scalar verify element 0
 *     3 = vle8.v masked: active+inactive elements
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_data
    vle8.v v8, (t1)
    SAVE_CSRS
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 4
    CHECK_CSRS_UNCHANGED

    SET_TEST_NUM 2
    la t1, result_buf
    lbu t2, 0(t1)
    li t3, 1
    FAIL_IF_NE t2, t3

    /* Masked load: mask=0b0101, tu/mu policy */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_vd
    vle8.v v8, (t1)
    la t1, tc3_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    la t1, tc3_mem
    vle8.v v8, (t1), v0.t
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_data:
    .byte 0x01, 0x02, 0x03, 0x04
.align 1
tc3_mask:
    .byte 5
.align 1
tc3_vd:
    .byte 0xdd, 0xdd, 0xdd, 0xdd
tc3_mem:
    .byte 0xaa, 0xbb, 0xcc, 0xdd
tc3_exp:
    .byte 0xaa, 0xdd, 0xcc, 0xdd

.align 4
result_buf:  .space 256
witness_buf: .space 256

