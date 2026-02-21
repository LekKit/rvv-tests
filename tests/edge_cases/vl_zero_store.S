/* Auto-generated test for vse32.v
 * vl=0 no-op: vse32.v must not write memory
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vl=0 store: memory modified
 *     2 = vl=0 store: vstart not zero
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vl=0 test: vse32.v — memory must be unchanged */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    li t1, 0xFFFFFFFF
    vmv.v.x v8, t1
    li t0, 0
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_target
    vse32.v v8, (t1)
    SET_TEST_NUM 2
    csrr t2, vstart
    FAIL_IF_NZ t2
    SET_TEST_NUM 1
    CHECK_MEM tc1_target, tc1_original, 16

    PASS_TEST

.data
.align 2
tc1_target:
    .word 0xcafebabe, 0xcafebabe, 0xcafebabe, 0xcafebabe
tc1_original:
    .word 0xcafebabe, 0xcafebabe, 0xcafebabe, 0xcafebabe

.align 4
result_buf:  .space 256
witness_buf: .space 256

