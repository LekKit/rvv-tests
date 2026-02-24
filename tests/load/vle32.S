/* Auto-generated test for vle32.v
 * Unit-stride load of 32-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vle32.v: basic load-store round-trip
 *     2 = vle32.v: scalar verify element 0
 *     3 = vle32.v masked: active+inactive elements
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_data
    vle32.v v8, (t1)
    SAVE_CSRS
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_data, 16
    CHECK_CSRS_UNCHANGED

    SET_TEST_NUM 2
    la t1, result_buf
    lwu t2, 0(t1)
    li t3, 1
    FAIL_IF_NE t2, t3

    /* Masked load: mask=0b0101, tu/mu policy */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_vd
    vle32.v v8, (t1)
    la t1, tc3_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    la t1, tc3_mem
    vle32.v v8, (t1), v0.t
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_data:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
.align 1
tc3_mask:
    .byte 5
.align 2
tc3_vd:
    .word 0x000000dd, 0x000000dd, 0x000000dd, 0x000000dd
tc3_mem:
    .word 0x000000aa, 0x000000bb, 0x000000cc, 0x000000dd
tc3_exp:
    .word 0x000000aa, 0x000000dd, 0x000000cc, 0x000000dd

.align 4
result_buf:  .space 256
witness_buf: .space 256

