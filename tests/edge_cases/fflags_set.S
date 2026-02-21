/* Auto-generated test for vfadd.vv
 * FP flags: verify fflags NX and OF bits set correctly
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = fflags: NX bit not set after inexact add
 *     2 = fflags: OF bit not set after overflow add
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Sub-test A: inexact — 1.0f + 1.5*2^-24 sets NX */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrw fflags, x0
    li t1, 0x3f800000
    vmv.v.x v16, t1
    li t1, 0x33c00000
    vmv.v.x v20, t1
    vfadd.vv v8, v16, v20
    SET_TEST_NUM 1
    csrr t0, fflags
    andi t0, t0, 1
    FAIL_IF_Z t0

    /* Sub-test B: overflow — FLT_MAX + FLT_MAX sets OF */
    csrw fflags, x0
    li t1, 0x7f7fffff
    vmv.v.x v16, t1
    vmv.v.x v20, t1
    vfadd.vv v8, v16, v20
    SET_TEST_NUM 2
    csrr t0, fflags
    andi t0, t0, 4
    FAIL_IF_Z t0

    PASS_TEST

.data

.align 4
result_buf:  .space 256
witness_buf: .space 256

