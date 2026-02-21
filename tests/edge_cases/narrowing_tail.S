/* Auto-generated test for vnsrl.wi
 * Narrowing tail undisturbed: vnsrl.wi e64→e32 vl=2 — verify active elements and tail preserved
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = narrowing-tail: active elements wrong
 *     2 = narrowing-tail: tail bytes modified
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Narrowing tail test: vnsrl.wi e64→e32, vl=2, tu/mu */
    vsetvli t0, x0, e32, m1, tu, mu
    li t1, 0xbeefdead
    vmv.v.x v8, t1
    li t0, 2
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_src64
    vle64.v v16, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    vnsrl.wi v8, v16, 0
    la t1, tc1_resbuf
    vs1r.v v8, (t1)
    SET_TEST_NUM 1
    CHECK_MEM tc1_resbuf, tc1_exp, 8
    SET_TEST_NUM 2
    la a1, tc1_resbuf
    addi a1, a1, 8
    la a2, tc1_bg
    addi a2, a2, 8
    csrr a3, vlenb
    addi a3, a3, -8
    blez a3, 1f
    jal ra, _mem_compare
    FAIL_IF_NZ a0
1:

    PASS_TEST

.data
.align 3
tc1_src64:
    .dword 0x0000000100000064, 0x00000002000000c8
.align 2
tc1_exp:
    .word 0x00000064, 0x000000c8
.align 4
tc1_bg:
    .fill 128, 4, 0xbeefdead
.align 4
tc1_resbuf:
    .space 512

.align 4
result_buf:  .space 256
witness_buf: .space 256

