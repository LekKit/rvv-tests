/* Auto-generated test for vwadd.vv
 * Full-VLMAX tail undisturbed: vwadd.vv e32->e64 — verifies ALL tail bytes in dest group
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlmax-tail widening: active elements wrong
 *     2 = vlmax-tail widening: tail elements modified
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Full-VLMAX tail test: vwadd.vv e32→e64, vl=2, tu/mu */
    /* Dest is v8-v9 (EMUL=2). Uses vs2r.v to dump, 2*vlenb for tail. */
    vsetvli t0, x0, e32, m2, tu, mu
    li t1, 0xdeadbeef
    vmv.v.x v8, t1
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    vwadd.vv v8, v16, v20
    la t1, tc1_resbuf
    vs2r.v v8, (t1)
    SET_TEST_NUM 1
    CHECK_MEM tc1_resbuf, tc1_exp, 16
    SET_TEST_NUM 2
    la a1, tc1_resbuf
    addi a1, a1, 16
    la a2, tc1_bg
    addi a2, a2, 16
    csrr a3, vlenb
    slli a3, a3, 1
    addi a3, a3, -16
    blez a3, 1f
    jal ra, _mem_compare
    FAIL_IF_NZ a0
1:

    PASS_TEST

.data
.align 2
tc1_s2:
    .word 0x00000064, 0x000000c8, 0x00000000, 0x00000000
tc1_s1:
    .word 0x0000000a, 0x00000014, 0x00000000, 0x00000000
.align 3
tc1_exp:
    .dword 0x000000000000006e, 0x00000000000000dc
.align 4
tc1_bg:
    .fill 256, 4, 0xdeadbeef
.align 4
tc1_resbuf:
    .space 1024

.align 4
result_buf:  .space 256
witness_buf: .space 256

