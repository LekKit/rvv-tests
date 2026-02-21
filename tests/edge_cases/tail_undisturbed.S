/* Auto-generated test for vadd.vv
 * Tail undisturbed: elements >= vl must be preserved
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = tail undisturbed: wrong result
 *     2 = tail undisturbed: vstart not zero after
 */
#include "riscv_test.h"
#include "test_macros.h"

    /* Test 1: vadd.vv with vl=2, tu — tail elements preserved */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_vd
    vle32.v v8, (t1)
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vadd.vv v8, v16, v20
    SET_TEST_NUM 2
    CHECK_VSTART_ZERO
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16

    PASS_TEST

.data
.align 2
tc1_vd:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
tc1_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
tc1_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc1_exp:
    .word 0x0000006e, 0x000000dc, 0x12345678, 0x9abcdef0

.align 4
result_buf:  .space 256
witness_buf: .space 256

