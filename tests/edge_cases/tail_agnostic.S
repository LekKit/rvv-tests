/* Auto-generated test for vadd.vv
 * Tail agnostic (ta): vl=2, verify only active elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = tail agnostic: active elements wrong
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Tail agnostic: vadd.vv e32 vl=2, ta/mu — only check active */
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, ta, mu
    vadd.vv v8, v16, v20
    la t1, tc1_resbuf
    vs1r.v v8, (t1)
    SET_TEST_NUM 1
    CHECK_MEM tc1_resbuf, tc1_exp, 8

    PASS_TEST

.data
.align 2
tc1_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
tc1_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc1_exp:
    .word 0x0000006e, 0x000000dc
.align 4
tc1_resbuf:
    .space 512

.align 4
result_buf:  .space 256
witness_buf: .space 256

