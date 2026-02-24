/* Auto-generated test for vsetvli/vsetivli edge cases
 * Tests AVL>VLMAX, rs1=x0 semantics, vsetivli range, and SEW/LMUL validity
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = AVL>VLMAX: vl not clipped to VLMAX
 *     2 = rs1=x0, rd!=x0: vl != VLMAX
 *     3 = rs1=x0, rd=x0: vl changed unexpectedly
 *     4 = vsetivli AVL=0: vl not 0
 *     5 = vsetivli AVL=31: vl wrong
 *     6 = valid SEW/LMUL combo set vill
 *     7 = invalid SEW/LMUL did NOT set vill
 *     8 = vsetvl rd, rs1, rs2: vl wrong
 *     9 = vsetvl rd, rs1, rs2: vadd with configured vtype wrong
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* AVL > VLMAX: vl should be clipped */
    li t0, 1000
    vsetvli t0, t0, e32, m1, ta, ma
    SET_TEST_NUM 1
    vsetvli t1, x0, e32, m1, ta, ma
    FAIL_IF_NE t0, t1

    /* rs1=x0, rd!=x0 → set vl=VLMAX */
    vsetvli t0, x0, e32, m1, ta, ma
    vsetvli t1, x0, e32, m1, ta, ma
    SET_TEST_NUM 2
    FAIL_IF_NE t0, t1

    /* rs1=x0, rd=x0 → keep vl (same SEW/LMUL ratio) */
    li t0, 3
    vsetvli t0, t0, e32, m1, ta, ma
    vsetvli x0, x0, e16, mf2, ta, ma
    csrr t1, vl
    SET_TEST_NUM 3
    li t2, 3
    FAIL_IF_NE t1, t2

    /* vsetivli AVL=0 */
    vsetivli t0, 0, e32, m1, ta, ma
    SET_TEST_NUM 4
    FAIL_IF_NZ t0

    /* vsetivli AVL=31 → vl = min(31, VLMAX) */
    vsetivli t0, 31, e32, m1, ta, ma
    vsetvli t1, x0, e32, m1, ta, ma
    SET_TEST_NUM 5
    li t2, 31
    bge t1, t2, 1f
    mv t2, t1
1:
    FAIL_IF_NE t0, t2

    /* All valid SEW/LMUL combos must not set vill */
    SET_TEST_NUM 6
    vsetvli t0, x0, e8, m1, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e8, m2, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e8, m4, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e8, m8, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e8, mf2, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e8, mf4, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e8, mf8, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e16, m1, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e16, m2, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e16, m4, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e16, m8, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e16, mf2, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e16, mf4, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e32, m1, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e32, m2, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e32, m4, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e32, m8, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e32, mf2, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e64, m1, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e64, m2, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e64, m4, ta, ma
    FAIL_IF_Z t0
    vsetvli t0, x0, e64, m8, ta, ma
    FAIL_IF_Z t0

    /* Invalid SEW/LMUL combos must set vill */
    SET_TEST_NUM 7
    /*   e32/mf4 (vtype=0xd6) */
    li t1, 214
    li t0, 4
    vsetvl t0, t0, t1
    FAIL_IF_NZ t0
    /*   e64/mf2 (vtype=0xdf) */
    li t1, 223
    li t0, 4
    vsetvl t0, t0, t1
    FAIL_IF_NZ t0
    /*   e16/mf8 (vtype=0xcd) */
    li t1, 205
    li t0, 4
    vsetvl t0, t0, t1
    FAIL_IF_NZ t0
    /*   e64/mf4 (vtype=0xde) */
    li t1, 222
    li t0, 4
    vsetvl t0, t0, t1
    FAIL_IF_NZ t0
    /*   e64/mf8 (vtype=0xdd) */
    li t1, 221
    li t0, 4
    vsetvl t0, t0, t1
    FAIL_IF_NZ t0
    /*   e32/mf8 (vtype=0xd5) */
    li t1, 213
    li t0, 4
    vsetvl t0, t0, t1
    FAIL_IF_NZ t0

    /* vsetvl register-register form */
    li t1, 208
    li t0, 4
    vsetvl t0, t0, t1
    SET_TEST_NUM 8
    li t2, 4
    FAIL_IF_NE t0, t2
    SET_TEST_NUM 9
    vmv.v.i v16, 10
    vmv.v.i v20, 5
    vadd.vv v8, v16, v20
    la t1, result_buf
    vse32.v v8, (t1)
    li t2, 15
    lw t3, 0(t1)
    FAIL_IF_NE t3, t2

    PASS_TEST

.data

.align 4
result_buf:  .space 256
witness_buf: .space 256

