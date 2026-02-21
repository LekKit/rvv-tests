/* Auto-generated test for vsadd.vv / vadd.vv
 * vxsat sticky: set by saturating op, preserved by non-saturating op, cleared only by explicit csrw
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vxsat: saturating result wrong
 *     2 = vxsat: not set after saturating add
 *     3 = vxsat: cleared by non-saturating add (not sticky)
 *     4 = vxsat: not cleared after csrw vxsat,0
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vxsat sticky test */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    csrw vxsat, x0
    li t1, 0x7FFFFFFF
    vmv.v.x v16, t1
    li t1, 1
    vmv.v.x v20, t1
    vsadd.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, vxsat_sat_exp, 16
    SET_TEST_NUM 2
    csrr t0, vxsat
    andi t0, t0, 1
    FAIL_IF_Z t0
    li t1, 1
    vmv.v.x v16, t1
    li t1, 2
    vmv.v.x v20, t1
    vadd.vv v8, v16, v20
    SET_TEST_NUM 3
    csrr t0, vxsat
    andi t0, t0, 1
    FAIL_IF_Z t0
    csrw vxsat, x0
    SET_TEST_NUM 4
    csrr t0, vxsat
    FAIL_IF_NZ t0

    PASS_TEST

.data
.align 2
vxsat_sat_exp:
    .word 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff

.align 4
result_buf:  .space 256
witness_buf: .space 256

