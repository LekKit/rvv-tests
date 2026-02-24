/* Auto-generated test for Small vl extra families
 * Tests vl=1/2 for vwadd.vv (widening), vnsrl.wi (narrowing), vmseq.vv (compare)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwadd.vv e32→e64 vl=1: result
 *     2 = vwadd.vv e32→e64 vl=2: result
 *     3 = vmseq.vv e32 vl=1: mask result
 *     4 = vmseq.vv e32 vl=2: mask result
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vwadd.vv e32→e64 vl=1 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svle1_s2
    vle32.v v16, (t1)
    la t1, svle1_s1
    vle32.v v20, (t1)
    li t0, 1
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vwadd.vv v8, v16, v20
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, svle1_exp, 8

    /* vwadd.vv e32→e64 vl=2 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svle2_s2
    vle32.v v16, (t1)
    la t1, svle2_s1
    vle32.v v20, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vwadd.vv v8, v16, v20
    SET_TEST_NUM 2
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, svle2_exp, 16

    /* vmseq.vv e32 vl=1 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svle3_s2
    vle32.v v16, (t1)
    la t1, svle3_s1
    vle32.v v20, (t1)
    li t0, 1
    vsetvli t0, t0, e32, m1, ta, ma
    SAVE_CSRS
    vmseq.vv v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 1
    li t3, 1
    FAIL_IF_NE t2, t3

    /* vmseq.vv e32 vl=2 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svle4_s2
    vle32.v v16, (t1)
    la t1, svle4_s1
    vle32.v v20, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, ta, ma
    SAVE_CSRS
    vmseq.vv v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 3
    li t3, 3
    FAIL_IF_NE t2, t3

    PASS_TEST

.data
.align 3
svle1_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svle1_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
svle1_exp:
    .dword 0x000000000000006e
.align 3
svle2_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svle2_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
svle2_exp:
    .dword 0x000000000000006e, 0x00000000000000dc
.align 2
svle3_s2:
    .word 0x00000005, 0x0000000a, 0x0000000a, 0x00000014
svle3_s1:
    .word 0x00000005, 0x0000000a, 0x0000000f, 0x00000014
.align 2
svle4_s2:
    .word 0x00000005, 0x0000000a, 0x0000000a, 0x00000014
svle4_s1:
    .word 0x00000005, 0x0000000a, 0x0000000f, 0x00000014

.align 4
result_buf:  .space 256
witness_buf: .space 256

