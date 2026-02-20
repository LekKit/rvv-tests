/* Auto-generated test for vsetvli
 * Set vector length and type from register + immediate
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsetvli e32/m1 avl=4: wrong vl
 *     2 = vsetvli: vl CSR mismatch
 *     3 = vsetvli: vtype.vill set
 *     4 = vsetvli e8/m1: vl>0
 *     5 = vsetvli e16/m1: vl>0
 *     6 = vsetvli e32/m1: vl>0
 *     7 = vsetvli e64/m1: vl>0
 *     8 = vsetvli e32/m1: no vill
 *     9 = vsetvli e32/m2: no vill
 *    10 = vsetvli e32/m4: no vill
 *    11 = vsetvli e32/m8: no vill
 *    12 = vsetvli e32/mf2: no vill
 *    13 = vsetvli e32/mf4: no vill
 *    14 = vsetvli e32/mf8: no vill
 *    15 = vsetvli avl=0: vl should be 0
 *    16 = vsetvli: vstart should be 0
 *    17 = vsetivli uimm=4 e32/m1: vl=4
 *    18 = vsetivli uimm=0: vl=0
 *    19 = vsetvl: basic
 *    20 = vsetvli: query vlmax
 *    21 = vsetvli: ta/ma encoding
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vsetvli e32, m1, tu, mu with avl=4 */
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t1, t0, e32, m1, tu, mu
    li t2, 4
    FAIL_IF_NE t1, t2
    SET_TEST_NUM 2
    csrr t3, vl
    FAIL_IF_NE t1, t3
    SET_TEST_NUM 3
    csrr t3, vtype
    srai t4, t3, 63
    FAIL_IF_NZ t4

    SET_TEST_NUM 4
    li t0, 4
    vsetvli t1, t0, e8, m1, tu, mu
    FAIL_IF_Z t1

    SET_TEST_NUM 5
    li t0, 4
    vsetvli t1, t0, e16, m1, tu, mu
    FAIL_IF_Z t1

    SET_TEST_NUM 6
    li t0, 4
    vsetvli t1, t0, e32, m1, tu, mu
    FAIL_IF_Z t1

    SET_TEST_NUM 7
    li t0, 4
    vsetvli t1, t0, e64, m1, tu, mu
    FAIL_IF_Z t1

    SET_TEST_NUM 8
    li t0, 4
    vsetvli t1, t0, e32, m1, tu, mu
    csrr t3, vtype

    SET_TEST_NUM 9
    li t0, 4
    vsetvli t1, t0, e32, m2, tu, mu
    csrr t3, vtype

    SET_TEST_NUM 10
    li t0, 4
    vsetvli t1, t0, e32, m4, tu, mu
    csrr t3, vtype

    SET_TEST_NUM 11
    li t0, 4
    vsetvli t1, t0, e32, m8, tu, mu
    csrr t3, vtype

    SET_TEST_NUM 12
    li t0, 4
    vsetvli t1, t0, e32, mf2, tu, mu
    csrr t3, vtype

    SET_TEST_NUM 13
    li t0, 4
    vsetvli t1, t0, e32, mf4, tu, mu
    csrr t3, vtype

    SET_TEST_NUM 14
    li t0, 4
    vsetvli t1, t0, e32, mf8, tu, mu
    csrr t3, vtype

    SET_TEST_NUM 15
    li t0, 0
    vsetvli t1, t0, e32, m1, tu, mu
    FAIL_IF_NZ t1

    SET_TEST_NUM 16
    li t0, 4
    vsetvli t1, t0, e32, m1, tu, mu
    csrr t3, vstart
    FAIL_IF_NZ t3

    SET_TEST_NUM 17
    vsetivli t1, 4, e32, m1, tu, mu
    li t2, 4
    FAIL_IF_NE t1, t2

    SET_TEST_NUM 18
    vsetivli t1, 0, e32, m1, tu, mu
    FAIL_IF_NZ t1

    SET_TEST_NUM 19
    li t0, 4
    vsetvli zero, t0, e32, m1, tu, mu
    csrr t2, vtype
    vsetvl t1, t0, t2
    li t3, 4
    FAIL_IF_NE t1, t3

    SET_TEST_NUM 20
    li t0, -1
    vsetvli t1, t0, e8, m8, tu, mu
    csrr t2, vlenb
    slli t2, t2, 3
    FAIL_IF_NE t1, t2

    SET_TEST_NUM 21
    li t0, 4
    vsetvli t1, t0, e32, m1, ta, ma
    csrr t3, vtype
    andi t4, t3, 0xC0
    li t5, 0xC0
    FAIL_IF_NE t4, t5

    PASS_TEST

.data

.align 4
result_buf:  .space 256
witness_buf: .space 256

