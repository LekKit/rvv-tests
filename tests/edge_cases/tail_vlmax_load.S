/* Auto-generated test for vle32.v
 * Full-VLMAX tail undisturbed: vle32.v — verifies ALL tail bytes up to VLEN/8
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlmax-tail load: active elements wrong
 *     2 = vlmax-tail load: tail elements modified
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Full-VLMAX tail test: vle32.v e32, vl=2, tu/mu */
    vsetvli t0, x0, e32, m1, tu, mu
    li t1, 0xdeadbeef
    vmv.v.x v8, t1
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    vle32.v v8, (t1)
    la t1, tc1_resbuf
    vs1r.v v8, (t1)
    SET_TEST_NUM 1
    CHECK_MEM tc1_resbuf, tc1_src, 8
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
tc1_src:
    .word 0x11111111, 0x22222222, 0x33333333, 0x44444444
.align 4
tc1_bg:
    .fill 128, 4, 0xdeadbeef
.align 4
tc1_resbuf:
    .space 512

.align 4
result_buf:  .space 256
witness_buf: .space 256

