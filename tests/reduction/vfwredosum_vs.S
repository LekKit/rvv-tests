/* Auto-generated test for vfwredosum.vs
 * Widening FP reduction vfwredosum.vs
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfwredosum.vs: result
 *     2 = vfwredosum.vs: CSR side-effect
 *     3 = vfwredosum.vs masked: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_vec
    vle32.v v16, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_init
    vle64.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vfwredosum.vs v8, v16, v20
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED_FP
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8

    /* Masked vfwredosum.vs */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_vec
    vle32.v v16, (t1)
    la t1, tc3_mask
    vlm.v v0, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc3_init
    vle64.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vfwredosum.vs v8, v16, v20, v0.t
    SET_TEST_NUM 3
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8

    PASS_TEST

.data
.align 3
tc1_vec:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc1_exp:
    .dword 0x4024000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 1
tc3_mask:
    .byte 5
.align 3
tc3_vec:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc3_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc3_exp:
    .dword 0x4010000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

