/* Auto-generated test for RVV 1.0 detection
 * Verifies the implementation is standard RVV 1.0 (not XTheadVector) by testing instructions that exist only in RVV 1.0: vsetivli, vlenb CSR, vmv1r.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsetivli not supported (XTheadVector?)
 *     2 = csrr vlenb failed
 *     3 = vlenb is zero (invalid)
 *     4 = vmv1r.v not supported (XTheadVector?)
 *     5 = vmv1r.v data integrity failed
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* RVV 1.0 detection: test instructions absent in XTheadVector */
    SET_TEST_NUM 1
    vsetivli t0, 4, e32, m1, ta, ma
    /* If we reach here, vsetivli works → not XTheadVector */
    SET_TEST_NUM 2
    csrr t0, vlenb
    SET_TEST_NUM 3
    FAIL_IF_Z t0
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, rvvdet_src
    vle32.v v16, (t1)
    SET_TEST_NUM 4
    vmv1r.v v8, v16
    SET_TEST_NUM 5
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, rvvdet_src, 16

    PASS_TEST

.data
.align 2
rvvdet_src:
    .word 0x11111111, 0x22222222, 0x33333333, 0x44444444

.align 4
result_buf:  .space 256
witness_buf: .space 256

