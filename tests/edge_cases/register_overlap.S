/* Auto-generated test for Register overlap constraints
 * Tests allowed register overlaps: vd==vs2, vd==vs1 for arithmetic; widening vd overlaps vs2; narrowing vs2 overlaps vd.
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vadd.vv vd==vs2: result wrong
 *     2 = vadd.vv vd==vs1: result wrong
 *     3 = vmul.vv vd==vs2==vs1 (square): result wrong
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vadd.vv v8, v8, v20 (vd == vs2) */
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    la t1, rovlp_s2
    vle32.v v8, (t1)
    la t1, rovlp_s1
    vle32.v v20, (t1)
    vadd.vv v8, v8, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, rovlp_exp_add, 16

    /* vadd.vv v8, v16, v8 (vd == vs1) */
    la t1, rovlp_s1
    vle32.v v8, (t1)
    la t1, rovlp_s2
    vle32.v v16, (t1)
    vadd.vv v8, v16, v8
    SET_TEST_NUM 2
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, rovlp_exp_add, 16

    /* vmul.vv v8, v8, v8 (square, all same register) */
    la t1, rovlp_sq
    vle32.v v8, (t1)
    vmul.vv v8, v8, v8
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, rovlp_exp_sq, 16

    PASS_TEST

.data
.align 4
rovlp_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
rovlp_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
rovlp_exp_add:
    .word 0x0000000b, 0x00000016, 0x00000021, 0x0000002c
rovlp_sq:
    .word 0x00000003, 0x00000005, 0x00000007, 0x0000000b
rovlp_exp_sq:
    .word 0x00000009, 0x00000019, 0x00000031, 0x00000079

.align 4
result_buf:  .space 256
witness_buf: .space 256

