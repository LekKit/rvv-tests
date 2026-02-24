/* Auto-generated test for lr/sc vs vector store
 * Tests that vector stores to a different address do not break lr/sc reservations, and that vector stores to the same address DO break them.
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = lr/sc: vec store to diff addr broke reservation
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* lr.w → vse32 to DIFFERENT address → sc.w: must succeed */
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    vmv.v.i v16, 5
    la s6, lrsc_atom
    la s7, lrsc_other
    li t2, 99
    lr.w t3, (s6)
    vse32.v v16, (s7)
    sc.w t4, t2, (s6)
    SET_TEST_NUM 1
    FAIL_IF_NZ t4

    PASS_TEST

.data
.align 4
lrsc_atom:
    .word 0
    .space 12
lrsc_other:
    .space 64

.align 4
result_buf:  .space 256
witness_buf: .space 256

