/* Auto-generated test for vfadd.vv
 * vl=0 no-op: vfadd.vv must not modify vd
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vl=0 fp: vd modified
 *     2 = vl=0 fp: vstart not zero
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vl=0 test: vfadd.vv — vd must be unchanged */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_bg
    vle32.v v8, (t1)
    la t1, tc1_src
    vle32.v v16, (t1)
    vle32.v v20, (t1)
    li t0, 0
    vsetvli t0, t0, e32, m1, tu, mu
    vfadd.vv v8, v16, v20
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
    .word 0x3f800000, 0x3f800000, 0x3f800000, 0x3f800000

.align 4
result_buf:  .space 256
witness_buf: .space 256

