/* Auto-generated test for vwredsumu.vs
 * Widening reduction vwredsumu.vs
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwredsumu.vs e8 basic: result
 *     2 = vwredsumu.vs e8 init: result
 *     3 = vwredsumu.vs e8 max: result
 *     4 = vwredsumu.vs e16 basic: result
 *     5 = vwredsumu.vs e16 init: result
 *     6 = vwredsumu.vs e16 max: result
 *     7 = vwredsumu.vs e32 basic: result
 *     8 = vwredsumu.vs e32 init: result
 *     9 = vwredsumu.vs e32 max: result
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
    SET_TEST_NUM 1
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 2
    CHECK_VSTART_ZERO

    /* Test 2: vwredsumu.vs SEW=8 init */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_vec
    vle8.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_init
    vle16.v v20, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 2
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 2
    CHECK_VSTART_ZERO

    /* Test 3: vwredsumu.vs SEW=8 max */
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
    SET_TEST_NUM 3
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 2
    CHECK_VSTART_ZERO

    /* Test 4: vwredsumu.vs SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc4_vec
    vle16.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc4_init
    vle32.v v20, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_VSTART_ZERO

    /* Test 5: vwredsumu.vs SEW=16 init */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_vec
    vle16.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc5_init
    vle32.v v20, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 5
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 4
    CHECK_VSTART_ZERO

    /* Test 6: vwredsumu.vs SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_vec
    vle16.v v16, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc6_init
    vle32.v v20, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 6
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 4
    CHECK_VSTART_ZERO

    /* Test 7: vwredsumu.vs SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc7_vec
    vle32.v v16, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc7_init
    vle64.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 7
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8
    CHECK_VSTART_ZERO

    /* Test 8: vwredsumu.vs SEW=32 init */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc8_vec
    vle32.v v16, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc8_init
    vle64.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 8
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 8
    CHECK_VSTART_ZERO

    /* Test 9: vwredsumu.vs SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_vec
    vle32.v v16, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc9_init
    vle64.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vwredsumu.vs v8, v16, v20
    SET_TEST_NUM 9
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 8
    CHECK_VSTART_ZERO

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
tc2_vec:
    .byte 0x01, 0x02, 0x03, 0x04
tc2_init:
    .half 0x0064, 0x0000, 0x0000, 0x0000
tc2_exp:
    .half 0x006e, 0x0000, 0x0000, 0x0000
.align 1
tc3_vec:
    .byte 0xff, 0xff, 0xff, 0xff
tc3_init:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc3_exp:
    .half 0x03fc, 0x0000, 0x0000, 0x0000
.align 2
tc4_vec:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc4_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc4_exp:
    .word 0x0000000a, 0x00000000, 0x00000000, 0x00000000
.align 2
tc5_vec:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc5_init:
    .word 0x00000064, 0x00000000, 0x00000000, 0x00000000
tc5_exp:
    .word 0x0000006e, 0x00000000, 0x00000000, 0x00000000
.align 2
tc6_vec:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc6_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc6_exp:
    .word 0x0003fffc, 0x00000000, 0x00000000, 0x00000000
.align 3
tc7_vec:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc7_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc7_exp:
    .dword 0x000000000000000a, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc8_vec:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc8_init:
    .dword 0x0000000000000064, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc8_exp:
    .dword 0x000000000000006e, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 3
tc9_vec:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc9_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc9_exp:
    .dword 0x00000003fffffffc, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

