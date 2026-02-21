/* Auto-generated test for vadd.vv
 * Fractional LMUL: e32/mf2, e16/mf4, e8/mf8 — verify active elements computed, tail bytes preserved
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = fract mf2: active elements wrong
 *     2 = fract mf2: tail bytes modified
 *     3 = fract mf4: active elements wrong
 *     4 = fract mf4: tail bytes modified
 *     5 = fract mf8: active elements wrong
 *     6 = fract mf8: tail bytes modified
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Sub-test: vadd.vv at e32/mf2 */
    vsetvli t0, x0, e32, m1, tu, mu
    li t1, 0xdeadbeef
    vmv.v.x v8, t1
    vsetvli s6, x0, e32, mf2, tu, mu
    li t1, 100
    vmv.v.x v16, t1
    li t1, 10
    vmv.v.x v20, t1
    vadd.vv v8, v16, v20
    la t1, fract_mf2_resbuf
    vs1r.v v8, (t1)
    SET_TEST_NUM 1
    la a1, fract_mf2_resbuf
    la a2, fract_mf2_exp
    slli a3, s6, 2
    jal ra, _mem_compare
    FAIL_IF_NZ a0
    SET_TEST_NUM 2
    la a1, fract_mf2_resbuf
    slli t0, s6, 2
    add a1, a1, t0
    la a2, fract_mf2_bg
    add a2, a2, t0
    csrr a3, vlenb
    sub a3, a3, t0
    blez a3, 1f
    jal ra, _mem_compare
    FAIL_IF_NZ a0
1:

    /* Sub-test: vadd.vv at e16/mf4 */
    vsetvli t0, x0, e32, m1, tu, mu
    li t1, 0xdeadbeef
    vmv.v.x v8, t1
    vsetvli s6, x0, e16, mf4, tu, mu
    li t1, 100
    vmv.v.x v16, t1
    li t1, 10
    vmv.v.x v20, t1
    vadd.vv v8, v16, v20
    la t1, fract_mf4_resbuf
    vs1r.v v8, (t1)
    SET_TEST_NUM 3
    la a1, fract_mf4_resbuf
    la a2, fract_mf4_exp
    slli a3, s6, 1
    jal ra, _mem_compare
    FAIL_IF_NZ a0
    SET_TEST_NUM 4
    la a1, fract_mf4_resbuf
    slli t0, s6, 1
    add a1, a1, t0
    la a2, fract_mf4_bg
    add a2, a2, t0
    csrr a3, vlenb
    sub a3, a3, t0
    blez a3, 1f
    jal ra, _mem_compare
    FAIL_IF_NZ a0
1:

    /* Sub-test: vadd.vv at e8/mf8 */
    vsetvli t0, x0, e32, m1, tu, mu
    li t1, 0xdeadbeef
    vmv.v.x v8, t1
    vsetvli s6, x0, e8, mf8, tu, mu
    li t1, 100
    vmv.v.x v16, t1
    li t1, 10
    vmv.v.x v20, t1
    vadd.vv v8, v16, v20
    la t1, fract_mf8_resbuf
    vs1r.v v8, (t1)
    SET_TEST_NUM 5
    la a1, fract_mf8_resbuf
    la a2, fract_mf8_exp
    mv a3, s6
    jal ra, _mem_compare
    FAIL_IF_NZ a0
    SET_TEST_NUM 6
    la a1, fract_mf8_resbuf
    mv t0, s6
    add a1, a1, t0
    la a2, fract_mf8_bg
    add a2, a2, t0
    csrr a3, vlenb
    sub a3, a3, t0
    blez a3, 1f
    jal ra, _mem_compare
    FAIL_IF_NZ a0
1:

    PASS_TEST

.data
.align 4
fract_mf2_exp:
    .fill 128, 4, 0x0000006e
.align 4
fract_mf2_bg:
    .fill 128, 4, 0xdeadbeef
.align 4
fract_mf2_resbuf:
    .space 512
.align 2
fract_mf4_exp:
    .fill 256, 2, 0x006e
.align 4
fract_mf4_bg:
    .fill 128, 4, 0xdeadbeef
.align 4
fract_mf4_resbuf:
    .space 512
.align 1
fract_mf8_exp:
    .fill 512, 1, 0x6e
.align 4
fract_mf8_bg:
    .fill 128, 4, 0xdeadbeef
.align 4
fract_mf8_resbuf:
    .space 512

.align 4
result_buf:  .space 256
witness_buf: .space 256

