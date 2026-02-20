/* Auto-generated test for vredmaxu.vs
 * Reduction vredmaxu.vs
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vredmaxu.vs e8 basic: wrong result
 *     2 = vredmaxu.vs e8 basic: CSR side-effect
 *     3 = vredmaxu.vs e8 init: wrong result
 *     4 = vredmaxu.vs e8 init: CSR side-effect
 *     5 = vredmaxu.vs e8 max: wrong result
 *     6 = vredmaxu.vs e8 max: CSR side-effect
 *     7 = vredmaxu.vs e16 basic: wrong result
 *     8 = vredmaxu.vs e16 basic: CSR side-effect
 *     9 = vredmaxu.vs e16 init: wrong result
 *    10 = vredmaxu.vs e16 init: CSR side-effect
 *    11 = vredmaxu.vs e16 max: wrong result
 *    12 = vredmaxu.vs e16 max: CSR side-effect
 *    13 = vredmaxu.vs e32 basic: wrong result
 *    14 = vredmaxu.vs e32 basic: CSR side-effect
 *    15 = vredmaxu.vs e32 init: wrong result
 *    16 = vredmaxu.vs e32 init: CSR side-effect
 *    17 = vredmaxu.vs e32 max: wrong result
 *    18 = vredmaxu.vs e32 max: CSR side-effect
 *    19 = vredmaxu.vs e64 basic: wrong result
 *    20 = vredmaxu.vs e64 basic: CSR side-effect
 *    21 = vredmaxu.vs e64 init: wrong result
 *    22 = vredmaxu.vs e64 init: CSR side-effect
 *    23 = vredmaxu.vs e64 max: wrong result
 *    24 = vredmaxu.vs e64 max: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vredmaxu.vs SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_vec
    vle8.v v16, (t1)
    la t1, tc1_init
    vle8.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 1
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    /* Test 3-4: vredmaxu.vs SEW=8 init */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_vec
    vle8.v v16, (t1)
    la t1, tc3_init
    vle8.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 1
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED

    /* Test 5-6: vredmaxu.vs SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_vec
    vle8.v v16, (t1)
    la t1, tc5_init
    vle8.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 5
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 1
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-8: vredmaxu.vs SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_vec
    vle16.v v16, (t1)
    la t1, tc7_init
    vle16.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 2
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED

    /* Test 9-10: vredmaxu.vs SEW=16 init */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_vec
    vle16.v v16, (t1)
    la t1, tc9_init
    vle16.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 9
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 2
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED

    /* Test 11-12: vredmaxu.vs SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_vec
    vle16.v v16, (t1)
    la t1, tc11_init
    vle16.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 11
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 2
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13-14: vredmaxu.vs SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_vec
    vle32.v v16, (t1)
    la t1, tc13_init
    vle32.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 4
    SET_TEST_NUM 14
    CHECK_CSRS_UNCHANGED

    /* Test 15-16: vredmaxu.vs SEW=32 init */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc15_vec
    vle32.v v16, (t1)
    la t1, tc15_init
    vle32.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 15
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 4
    SET_TEST_NUM 16
    CHECK_CSRS_UNCHANGED

    /* Test 17-18: vredmaxu.vs SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc17_vec
    vle32.v v16, (t1)
    la t1, tc17_init
    vle32.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 17
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 4
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    /* Test 19-20: vredmaxu.vs SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc19_vec
    vle64.v v16, (t1)
    la t1, tc19_init
    vle64.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 19
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 8
    SET_TEST_NUM 20
    CHECK_CSRS_UNCHANGED

    /* Test 21-22: vredmaxu.vs SEW=64 init */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc21_vec
    vle64.v v16, (t1)
    la t1, tc21_init
    vle64.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 21
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc21_exp, 8
    SET_TEST_NUM 22
    CHECK_CSRS_UNCHANGED

    /* Test 23-24: vredmaxu.vs SEW=64 max */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc23_vec
    vle64.v v16, (t1)
    la t1, tc23_init
    vle64.v v20, (t1)
    SAVE_CSRS
    vredmaxu.vs v8, v16, v20
    SET_TEST_NUM 23
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc23_exp, 8
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_vec:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_init:
    .byte 0x00, 0x00, 0x00, 0x00
tc1_exp:
    .byte 0x04, 0x00, 0x00, 0x00
.align 1
tc3_vec:
    .byte 0x01, 0x02, 0x03, 0x04
tc3_init:
    .byte 0x64, 0x00, 0x00, 0x00
tc3_exp:
    .byte 0x64, 0x00, 0x00, 0x00
.align 1
tc5_vec:
    .byte 0xff, 0xfe, 0xfd, 0x00
tc5_init:
    .byte 0x00, 0x00, 0x00, 0x00
tc5_exp:
    .byte 0xff, 0x00, 0x00, 0x00
.align 1
tc7_vec:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc7_init:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc7_exp:
    .half 0x0004, 0x0000, 0x0000, 0x0000
.align 1
tc9_vec:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc9_init:
    .half 0x0064, 0x0000, 0x0000, 0x0000
tc9_exp:
    .half 0x0064, 0x0000, 0x0000, 0x0000
.align 1
tc11_vec:
    .half 0xffff, 0xfffe, 0xfffd, 0x0000
tc11_init:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc11_exp:
    .half 0xffff, 0x0000, 0x0000, 0x0000
.align 2
tc13_vec:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc13_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc13_exp:
    .word 0x00000004, 0x00000000, 0x00000000, 0x00000000
.align 2
tc15_vec:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc15_init:
    .word 0x00000064, 0x00000000, 0x00000000, 0x00000000
tc15_exp:
    .word 0x00000064, 0x00000000, 0x00000000, 0x00000000
.align 2
tc17_vec:
    .word 0xffffffff, 0xfffffffe, 0xfffffffd, 0x00000000
tc17_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc17_exp:
    .word 0xffffffff, 0x00000000, 0x00000000, 0x00000000
.align 3
tc19_vec:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc19_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc19_exp:
    .dword 0x0000000000000004, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc21_vec:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc21_init:
    .dword 0x0000000000000064, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc21_exp:
    .dword 0x0000000000000064, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc23_vec:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0xfffffffffffffffd, 0x0000000000000000
tc23_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc23_exp:
    .dword 0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

