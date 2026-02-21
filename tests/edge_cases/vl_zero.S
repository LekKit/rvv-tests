/* Auto-generated test for vadd.vv
 * vl=0 no-op: instruction must not modify vd
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vl=0: vd modified
 *     2 = vl=0: vstart not zero after
 */
#include "riscv_test.h"
#include "test_macros.h"

    /* Test 1: vadd.vv with vl=0 — vd must be unchanged */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_bg
    vle32.v v8, (t1)
    la t1, tc1_src
    vle32.v v16, (t1)
    vle32.v v20, (t1)
    li t0, 0
    vsetvli t0, t0, e32, m1, tu, mu
    vadd.vv v8, v16, v20
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
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004

.align 4
result_buf:  .space 256
witness_buf: .space 256

