/* Auto-generated test for vfwadd.wv
 * Widening FP vfwadd.wv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfwadd.wv basic: result
 *     2 = vfwadd.wv basic: CSR side-effect
 *     3 = vfwadd.wv neg: result
 *     4 = vfwadd.wv neg: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc1_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vfwadd.wv v8, v16, v20
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED_FP
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32

    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc3_s2w
    vle64.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vfwadd.wv v8, v16, v20
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED_FP
    SET_TEST_NUM 3
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 32

    PASS_TEST

.data
.align 3
tc1_s2w:
    .dword 0x4024000000000000, 0x4034000000000000, 0x403e000000000000, 0x4044000000000000
tc1_s1:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_exp:
    .dword 0x4026000000000000, 0x4036000000000000, 0x4040800000000000, 0x4046000000000000
.align 3
tc3_s2w:
    .dword 0xbff0000000000000, 0xc000000000000000, 0x4008000000000000, 0x4010000000000000
tc3_s1:
    .word 0x3f800000, 0x40000000, 0xc0400000, 0xc0800000
tc3_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

