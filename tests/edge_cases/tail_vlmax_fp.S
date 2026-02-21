/* Auto-generated test for vfadd.vv
 * Full-VLMAX tail undisturbed: vfadd.vv e32 — verifies ALL tail bytes up to VLEN/8
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlmax-tail fp: active elements wrong
 *     2 = vlmax-tail fp: tail elements modified
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Full-VLMAX tail test: vfadd.vv e32, vl=2, tu/mu */
    vsetvli t0, x0, e32, m1, tu, mu
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
    vfadd.vv v8, v16, v20
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
.align 2
tc1_s2:
    .word 0x3f800000, 0x40000000, 0x00000000, 0x00000000
tc1_s1:
    .word 0x3f000000, 0x3e800000, 0x00000000, 0x00000000
tc1_exp:
    .word 0x3fc00000, 0x40100000
.align 4
tc1_bg:
    .fill 128, 4, 0xdeadbeef
.align 4
tc1_resbuf:
    .space 512

.align 4
result_buf:  .space 256
witness_buf: .space 256

