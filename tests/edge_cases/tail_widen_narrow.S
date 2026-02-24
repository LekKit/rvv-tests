/* Auto-generated test for Tail undisturbed widening/narrowing
 * Verifies tail preservation (tu) for vwadd.vv (widening) and vnsrl.wi (narrowing) with vl=2
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwadd.vv e32→e64 tu vl=2: tail not preserved
 *     2 = vnsrl.wi e64→e32 tu vl=2: tail not preserved
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Widening tail: vwadd.vv e32→e64, vl=2, tu */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, twn1_pf64
    vle64.v v8, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, twn1_s2
    vle32.v v16, (t1)
    la t1, twn1_s1
    vle32.v v20, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vwadd.vv v8, v16, v20
    la t1, twn1_res
    vs2r.v v8, (t1)
    SET_TEST_NUM 1
    CHECK_MEM twn1_res, twn1_exp, 32

    /* Narrowing tail: vnsrl.wi e64→e32, vl=2, tu */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, twn2_pf32
    vle32.v v8, (t1)
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, twn2_s64
    vle64.v v16, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vnsrl.wi v8, v16, 0
    la t1, twn2_res
    vs1r.v v8, (t1)
    SET_TEST_NUM 2
    CHECK_MEM twn2_res, twn2_exp, 16

    PASS_TEST

.data
.align 4
twn1_pf64:
    .dword 0xaaaaaaaaaaaaaaaa, 0xbbbbbbbbbbbbbbbb, 0xcccccccccccccccc, 0xdddddddddddddddd
twn1_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
twn1_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
twn1_exp:
    .dword 0x000000000000006e, 0x00000000000000dc, 0xcccccccccccccccc, 0xdddddddddddddddd
twn1_res:
    .space 128
.align 4
twn2_pf32:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
twn2_s64:
    .dword 0x0000000100000002, 0x0000000300000004, 0x0000000500000006, 0x0000000700000008
twn2_exp:
    .word 0x00000002, 0x00000004, 0xcccccccc, 0xdddddddd
twn2_res:
    .space 64

.align 4
result_buf:  .space 256
witness_buf: .space 256

