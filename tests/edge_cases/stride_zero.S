/* Auto-generated test for vlse32.v
 * Stride=0 load: all elements read from same address (broadcast)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = stride-zero: wrong result
 *     2 = stride-zero: witness changed
 *     3 = stride-zero: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Strided load with stride=0 — broadcast single word */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_w
    vle32.v v24, (t1)
    SAVE_CSRS
    la t1, tc1_src
    vlse32.v v8, (t1), x0
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    SET_TEST_NUM 2
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 16
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0xcafebabe
tc1_exp:
    .word 0xcafebabe, 0xcafebabe, 0xcafebabe, 0xcafebabe
tc1_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0

.align 4
result_buf:  .space 256
witness_buf: .space 256

