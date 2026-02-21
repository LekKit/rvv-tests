/* Auto-generated test for vle32.v
 * vl=0 no-op: vle32.v must not modify vd
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vl=0 load: vd modified
 *     2 = vl=0 load: vstart not zero
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vl=0 test: vle32.v — vd must be unchanged */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_bg
    vle32.v v8, (t1)
    li t0, 0
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    vle32.v v8, (t1)
    SET_TEST_NUM 2
    csrr t2, vstart
    FAIL_IF_NZ t2
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_bg, 16

    PASS_TEST

.data
.align 2
tc1_bg:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
tc1_src:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff

.align 4
result_buf:  .space 256
witness_buf: .space 256

