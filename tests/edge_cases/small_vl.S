/* Auto-generated test for Small vl edge cases
 * Verifies correct behavior with vl=1, vl=2, vl=3 for representative instructions (vadd, vmul, vand, vle32, vse32, vslidedown, vredsum). Tail elements must be preserved (tu).
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vadd.vv vl=1: result
 *     2 = vmul.vv vl=1: result
 *     3 = vle32.v vl=1: tail preserved
 *     4 = vredsum.vs vl=1: result
 *     5 = vadd.vv vl=2: result
 *     6 = vmul.vv vl=2: result
 *     7 = vle32.v vl=2: tail preserved
 *     8 = vredsum.vs vl=2: result
 *     9 = vadd.vv vl=3: result
 *    10 = vmul.vv vl=3: result
 *    11 = vle32.v vl=3: tail preserved
 *    12 = vredsum.vs vl=3: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vadd.vv with vl=1 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl1_pf
    vle32.v v8, (t1)
    la t1, svl1_s2
    vle32.v v16, (t1)
    la t1, svl1_s1
    vle32.v v20, (t1)
    li t0, 1
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vadd.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl1_exp, 16

    /* vmul.vv with vl=1 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl2_pf
    vle32.v v8, (t1)
    la t1, svl2_s2
    vle32.v v16, (t1)
    la t1, svl2_s1
    vle32.v v20, (t1)
    li t0, 1
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vmul.vv v8, v16, v20
    SET_TEST_NUM 2
    la t1, result_buf
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl2_exp, 16

    /* vle32.v with vl=1 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl3_pf
    vle32.v v8, (t1)
    li t0, 1
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    la t1, svl3_mem
    vle32.v v8, (t1)
    SET_TEST_NUM 3
    la t1, result_buf
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl3_exp, 16

    /* vredsum.vs with vl=1 */
    li t0, 1
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl4_vec
    vle32.v v16, (t1)
    la t1, svl4_init
    vle32.v v20, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl4_exp, 4

    /* vadd.vv with vl=2 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl5_pf
    vle32.v v8, (t1)
    la t1, svl5_s2
    vle32.v v16, (t1)
    la t1, svl5_s1
    vle32.v v20, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vadd.vv v8, v16, v20
    SET_TEST_NUM 5
    la t1, result_buf
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl5_exp, 16

    /* vmul.vv with vl=2 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl6_pf
    vle32.v v8, (t1)
    la t1, svl6_s2
    vle32.v v16, (t1)
    la t1, svl6_s1
    vle32.v v20, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vmul.vv v8, v16, v20
    SET_TEST_NUM 6
    la t1, result_buf
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl6_exp, 16

    /* vle32.v with vl=2 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl7_pf
    vle32.v v8, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    la t1, svl7_mem
    vle32.v v8, (t1)
    SET_TEST_NUM 7
    la t1, result_buf
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl7_exp, 16

    /* vredsum.vs with vl=2 */
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl8_vec
    vle32.v v16, (t1)
    la t1, svl8_init
    vle32.v v20, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 8
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl8_exp, 4

    /* vadd.vv with vl=3 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl9_pf
    vle32.v v8, (t1)
    la t1, svl9_s2
    vle32.v v16, (t1)
    la t1, svl9_s1
    vle32.v v20, (t1)
    li t0, 3
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vadd.vv v8, v16, v20
    SET_TEST_NUM 9
    la t1, result_buf
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl9_exp, 16

    /* vmul.vv with vl=3 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl10_pf
    vle32.v v8, (t1)
    la t1, svl10_s2
    vle32.v v16, (t1)
    la t1, svl10_s1
    vle32.v v20, (t1)
    li t0, 3
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vmul.vv v8, v16, v20
    SET_TEST_NUM 10
    la t1, result_buf
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl10_exp, 16

    /* vle32.v with vl=3 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl11_pf
    vle32.v v8, (t1)
    li t0, 3
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    la t1, svl11_mem
    vle32.v v8, (t1)
    SET_TEST_NUM 11
    la t1, result_buf
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl11_exp, 16

    /* vredsum.vs with vl=3 */
    li t0, 3
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, svl12_vec
    vle32.v v16, (t1)
    la t1, svl12_init
    vle32.v v20, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 12
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, svl12_exp, 4

    PASS_TEST

.data
.align 2
svl1_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
svl1_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svl1_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
svl1_exp:
    .word 0x00000065, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
.align 2
svl2_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
svl2_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svl2_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
svl2_exp:
    .word 0x00000064, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
.align 2
svl3_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
svl3_mem:
    .word 0x00001111, 0x00002222, 0x00003333, 0x00004444
svl3_exp:
    .word 0x00001111, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
.align 2
svl4_vec:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svl4_init:
    .word 0x000003e8, 0x00000000, 0x00000000, 0x00000000
svl4_exp:
    .word 0x0000044c, 0x00000000, 0x00000000, 0x00000000
.align 2
svl5_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
svl5_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svl5_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
svl5_exp:
    .word 0x00000065, 0x000000ca, 0xcccccccc, 0xdddddddd
.align 2
svl6_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
svl6_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svl6_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
svl6_exp:
    .word 0x00000064, 0x00000190, 0xcccccccc, 0xdddddddd
.align 2
svl7_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
svl7_mem:
    .word 0x00001111, 0x00002222, 0x00003333, 0x00004444
svl7_exp:
    .word 0x00001111, 0x00002222, 0xcccccccc, 0xdddddddd
.align 2
svl8_vec:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svl8_init:
    .word 0x000003e8, 0x00000000, 0x00000000, 0x00000000
svl8_exp:
    .word 0x00000514, 0x00000000, 0x00000000, 0x00000000
.align 2
svl9_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
svl9_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svl9_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
svl9_exp:
    .word 0x00000065, 0x000000ca, 0x0000012f, 0xdddddddd
.align 2
svl10_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
svl10_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svl10_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
svl10_exp:
    .word 0x00000064, 0x00000190, 0x00000384, 0xdddddddd
.align 2
svl11_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
svl11_mem:
    .word 0x00001111, 0x00002222, 0x00003333, 0x00004444
svl11_exp:
    .word 0x00001111, 0x00002222, 0x00003333, 0xdddddddd
.align 2
svl12_vec:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
svl12_init:
    .word 0x000003e8, 0x00000000, 0x00000000, 0x00000000
svl12_exp:
    .word 0x00000640, 0x00000000, 0x00000000, 0x00000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

