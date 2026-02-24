/* Auto-generated test for vwredsumu.vs
 * Widening reduction vwredsumu.vs
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwredsumu.vs e8 basic: result
 *     2 = vwredsumu.vs e8 basic: CSR side-effect
 *     3 = vwredsumu.vs e8 init: result
 *     4 = vwredsumu.vs e8 init: CSR side-effect
 *     5 = vwredsumu.vs e8 max: result
 *     6 = vwredsumu.vs e8 max: CSR side-effect
 *     7 = vwredsumu.vs e16 basic: result
 *     8 = vwredsumu.vs e16 basic: CSR side-effect
 *     9 = vwredsumu.vs e16 init: result
 *    10 = vwredsumu.vs e16 init: CSR side-effect
 *    11 = vwredsumu.vs e16 max: result
 *    12 = vwredsumu.vs e16 max: CSR side-effect
 *    13 = vwredsumu.vs e32 basic: result
 *    14 = vwredsumu.vs e32 basic: CSR side-effect
 *    15 = vwredsumu.vs e32 init: result
 *    16 = vwredsumu.vs e32 init: CSR side-effect
 *    17 = vwredsumu.vs e32 max: result
 *    18 = vwredsumu.vs e32 max: CSR side-effect
 *    19 = vwredsumu.vs e8 masked: result
 *    20 = vwredsumu.vs e16 masked: result
 *    21 = vwredsumu.vs e32 masked: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vwredsumu.vs SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_vec
    vle8.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_init
    vle16.v v20, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 2

    /* Test 3: vwredsumu.vs SEW=8 init */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_vec
    vle8.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc3_init
    vle16.v v20, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 3
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 2

    /* Test 5: vwredsumu.vs SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_vec
    vle8.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_init
    vle16.v v20, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 5
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 2

    /* Test 7: vwredsumu.vs SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_vec
    vle16.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc7_init
    vle32.v v20, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 7
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 4

    /* Test 9: vwredsumu.vs SEW=16 init */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_vec
    vle16.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_init
    vle32.v v20, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 9
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 4

    /* Test 11: vwredsumu.vs SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_vec
    vle16.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_init
    vle32.v v20, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 11
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 4

    /* Test 13: vwredsumu.vs SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_vec
    vle32.v v16, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc13_init
    vle64.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 14
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 13
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 8

    /* Test 15: vwredsumu.vs SEW=32 init */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc15_vec
    vle32.v v16, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc15_init
    vle64.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 16
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 15
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 8

    /* Test 17: vwredsumu.vs SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc17_vec
    vle32.v v16, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc17_init
    vle64.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED
    SET_TEST_NUM 17
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 8

    /* Masked vwredsumu.vs SEW=8 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc19_vec
    vle8.v v16, (t1)
    la t1, tc19_mask
    vlm.v v0, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc19_init
    vle16.v v20, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20, v0.t
    SET_TEST_NUM 19
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 2

    /* Masked vwredsumu.vs SEW=16 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc20_vec
    vle16.v v16, (t1)
    la t1, tc20_mask
    vlm.v v0, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc20_init
    vle32.v v20, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20, v0.t
    SET_TEST_NUM 20
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc20_exp, 4

    /* Masked vwredsumu.vs SEW=32 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc21_vec
    vle32.v v16, (t1)
    la t1, tc21_mask
    vlm.v v0, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc21_init
    vle64.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20, v0.t
    SET_TEST_NUM 21
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc21_exp, 8

    PASS_TEST

.data
.align 1
tc1_vec:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_init:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc1_exp:
    .half 0x000a, 0x0000, 0x0000, 0x0000
.align 1
tc3_vec:
    .byte 0x01, 0x02, 0x03, 0x04
tc3_init:
    .half 0x0064, 0x0000, 0x0000, 0x0000
tc3_exp:
    .half 0x006e, 0x0000, 0x0000, 0x0000
.align 1
tc5_vec:
    .byte 0xff, 0xff, 0xff, 0xff
tc5_init:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc5_exp:
    .half 0x03fc, 0x0000, 0x0000, 0x0000
.align 2
tc7_vec:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc7_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc7_exp:
    .word 0x0000000a, 0x00000000, 0x00000000, 0x00000000
.align 2
tc9_vec:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc9_init:
    .word 0x00000064, 0x00000000, 0x00000000, 0x00000000
tc9_exp:
    .word 0x0000006e, 0x00000000, 0x00000000, 0x00000000
.align 2
tc11_vec:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc11_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc11_exp:
    .word 0x0003fffc, 0x00000000, 0x00000000, 0x00000000
.align 3
tc13_vec:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc13_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc13_exp:
    .dword 0x000000000000000a, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc15_vec:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc15_init:
    .dword 0x0000000000000064, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc15_exp:
    .dword 0x000000000000006e, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc17_vec:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc17_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc17_exp:
    .dword 0x00000003fffffffc, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 1
tc19_mask:
    .byte 5
.align 1
tc19_vec:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc19_init:
    .half 0x0064, 0x0000, 0x0000, 0x0000
tc19_exp:
    .half 0x008c, 0x0000, 0x0000, 0x0000
.align 1
tc20_mask:
    .byte 5
.align 2
tc20_vec:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc20_init:
    .word 0x00000064, 0x00000000, 0x00000000, 0x00000000
tc20_exp:
    .word 0x0000008c, 0x00000000, 0x00000000, 0x00000000
.align 1
tc21_mask:
    .byte 5
.align 3
tc21_vec:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc21_init:
    .dword 0x0000000000000064, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc21_exp:
    .dword 0x000000000000008c, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

