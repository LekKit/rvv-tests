/* Auto-generated test for vmv.x.s
 * Move vector element 0 to scalar register
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmv.x.s e8: result
 *     2 = vmv.x.s e16: result
 *     3 = vmv.x.s e32: result
 *     4 = vmv.x.s e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    li t2, 0x42
    vmv.v.x v16, t2
    SAVE_CSRS
    vmv.x.s t3, v16
    SET_TEST_NUM 1
    li t4, 0x42
    FAIL_IF_NE t3, t4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    li t2, 0x0042
    vmv.v.x v16, t2
    SAVE_CSRS
    vmv.x.s t3, v16
    SET_TEST_NUM 2
    li t4, 0x0042
    FAIL_IF_NE t3, t4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    li t2, 0x00000042
    vmv.v.x v16, t2
    SAVE_CSRS
    vmv.x.s t3, v16
    SET_TEST_NUM 3
    li t4, 0x00000042
    FAIL_IF_NE t3, t4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    li t2, 0x0000000000000042
    vmv.v.x v16, t2
    SAVE_CSRS
    vmv.x.s t3, v16
    SET_TEST_NUM 4
    li t4, 0x0000000000000042
    FAIL_IF_NE t3, t4
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data

.align 4
result_buf:  .space 256
witness_buf: .space 256

