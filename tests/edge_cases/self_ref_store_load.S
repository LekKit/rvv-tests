/* Auto-generated test for Self-referential store-load
 * Same vector register stores data, then loads it back from the same address. Tests that store retirement is correct.
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = self-ref v8 store-load: data mismatch
 *     2 = self-ref overwrite-load: got stale data
 *     3 = self-ref chain × 4: mismatch
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* v8 stores to buf, then v8 loads from same buf */
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    la t1, sref_data1
    vle32.v v8, (t1)
    la s6, sref_buf1
    vse32.v v8, (s6)
    vle32.v v8, (s6)
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, sref_data1, 16

    /* v8=A, store A, v8=B (via vadd), load from buf → must get A */
    la t1, sref_data_a
    vle32.v v8, (t1)
    la s6, sref_buf2
    vse32.v v8, (s6)
    la t1, sref_data_b
    vle32.v v16, (t1)
    vadd.vv v8, v16, v16
    vle32.v v8, (s6)
    SET_TEST_NUM 2
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, sref_data_a, 16

    /* Self-ref chain: 4 iterations, new data each time */
    la s6, sref_buf3
    SET_TEST_NUM 3
    li s7, 0
sref_loop:
    vmv.v.x v8, s7
    vid.v v16
    vadd.vv v8, v8, v16
    vse32.v v8, (s6)
    vle32.v v8, (s6)
    vmv.v.x v16, s7
    vid.v v20
    vadd.vv v16, v16, v20
    vmseq.vv v0, v8, v16
    vcpop.m t3, v0
    li t4, 4
    bne t3, t4, sref_fail
    addi s7, s7, 100
    li t4, 400
    blt s7, t4, sref_loop
    j sref_done
sref_fail:
    FAIL_TEST
sref_done:

    PASS_TEST

.data
.align 4
sref_data1:
    .word 0xaaaa1111, 0xbbbb2222, 0xcccc3333, 0xdddd4444
sref_data_a:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
sref_data_b:
    .word 0x000003e7, 0x000003e7, 0x000003e7, 0x000003e7
sref_buf1:
    .space 16
sref_buf2:
    .space 16
sref_buf3:
    .space 16

.align 4
result_buf:  .space 256
witness_buf: .space 256

