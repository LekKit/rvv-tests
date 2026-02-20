/* Auto-generated test for vmv.v.x
 * Scalar broadcast to vector
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmv.v.x e8: result
 *     2 = vmv.v.x e16: result
 *     3 = vmv.v.x e32: result
 *     4 = vmv.v.x e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    li a0, 0x2a
    SAVE_CSRS
    vmv.v.x v8, a0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    li a0, 0x002a
    SAVE_CSRS
    vmv.v.x v8, a0
    SET_TEST_NUM 2
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    li a0, 0x0000002a
    SAVE_CSRS
    vmv.v.x v8, a0
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    li a0, 0x000000000000002a
    SAVE_CSRS
    vmv.v.x v8, a0
    SET_TEST_NUM 4
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_exp:
    .byte 0x2a, 0x2a, 0x2a, 0x2a
.align 1
tc2_exp:
    .half 0x002a, 0x002a, 0x002a, 0x002a
.align 2
tc3_exp:
    .word 0x0000002a, 0x0000002a, 0x0000002a, 0x0000002a
.align 3
tc4_exp:
    .dword 0x000000000000002a, 0x000000000000002a, 0x000000000000002a, 0x000000000000002a

.align 4
result_buf:  .space 256
witness_buf: .space 256

