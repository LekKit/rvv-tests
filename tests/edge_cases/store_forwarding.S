/* Auto-generated test for Store-to-load forwarding
 * Tests that vector store-to-load forwarding works without fence (same hart, same address). RVWMO guarantees visibility.
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vse32→vle32 same address no fence: mismatch
 *     2 = vse32→vadd→vle32 same address: mismatch
 *     3 = back-to-back store-load × 4: mismatch
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Store then load, same address, no fence */
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    la t1, sfwd_data1
    vle32.v v16, (t1)
    la s6, sfwd_buf1
    vse32.v v16, (s6)
    vle32.v v8, (s6)
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, sfwd_data1, 16

    /* Store, intervening vadd, then load same address */
    la t1, sfwd_data2
    vle32.v v16, (t1)
    la s6, sfwd_buf2
    vse32.v v16, (s6)
    vadd.vv v20, v16, v16
    vle32.v v8, (s6)
    SET_TEST_NUM 2
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, sfwd_data2, 16

    /* 4 back-to-back store-load pairs, same address */
    la s6, sfwd_buf3
    SET_TEST_NUM 3
    li s7, 0
sfwd_loop:
    vmv.v.x v16, s7
    vid.v v20
    vadd.vv v16, v16, v20
    vse32.v v16, (s6)
    vle32.v v8, (s6)
    vmseq.vv v0, v8, v16
    vcpop.m t3, v0
    li t4, 4
    bne t3, t4, sfwd_fail
    addi s7, s7, 1
    li t4, 4
    blt s7, t4, sfwd_loop
    j sfwd_done
sfwd_fail:
    FAIL_TEST
sfwd_done:

    PASS_TEST

.data
.align 4
sfwd_data1:
    .word 0xdead0001, 0xdead0002, 0xdead0003, 0xdead0004
sfwd_data2:
    .word 0xcafe0001, 0xcafe0002, 0xcafe0003, 0xcafe0004
sfwd_buf1:
    .space 16
sfwd_buf2:
    .space 16
sfwd_buf3:
    .space 16

.align 4
result_buf:  .space 256
witness_buf: .space 256

