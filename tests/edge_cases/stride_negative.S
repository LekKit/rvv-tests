/* Auto-generated test for vlse32.v
 * Negative stride load: elements loaded in reverse order
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = stride-negative: wrong result
 *     2 = stride-negative: witness changed
 *     3 = stride-negative: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Strided load with stride=-4 — reverse order */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_w
    vle32.v v24, (t1)
    SAVE_CSRS
    la t1, tc1_src
    addi t1, t1, 12
    li a0, -4
    vlse32.v v8, (t1), a0
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
    .word 0xaaaa0001, 0xbbbb0002, 0xcccc0003, 0xdddd0004
tc1_exp:
    .word 0xdddd0004, 0xcccc0003, 0xbbbb0002, 0xaaaa0001
tc1_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0

.align 4
result_buf:  .space 256
witness_buf: .space 256

