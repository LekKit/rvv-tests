/* Auto-generated test for vnmsac.vx
 * Multiply-accumulate vnmsac.vx
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vnmsac.vx e8 basic: result
 *     2 = vnmsac.vx e8 accum: result
 *     3 = vnmsac.vx e16 basic: result
 *     4 = vnmsac.vx e16 accum: result
 *     5 = vnmsac.vx e32 basic: result
 *     6 = vnmsac.vx e32 accum: result
 *     7 = vnmsac.vx e64 basic: result
 *     8 = vnmsac.vx e64 accum: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vnmsac.vx SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_vd
    vle8.v v8, (t1)
    la t1, tc1_s2
    vle8.v v16, (t1)
    li a0, 0x02
    SAVE_CSRS
    vnmsac.vx v8, a0, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    /* Test 2: vnmsac.vx SEW=8 accum */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_vd
    vle8.v v8, (t1)
    la t1, tc2_s2
    vle8.v v16, (t1)
    li a0, 0x01
    SAVE_CSRS
    vnmsac.vx v8, a0, v16
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    /* Test 3: vnmsac.vx SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc3_vd
    vle16.v v8, (t1)
    la t1, tc3_s2
    vle16.v v16, (t1)
    li a0, 0x0002
    SAVE_CSRS
    vnmsac.vx v8, a0, v16
    SET_TEST_NUM 3
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    CHECK_CSRS_UNCHANGED

    /* Test 4: vnmsac.vx SEW=16 accum */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc4_vd
    vle16.v v8, (t1)
    la t1, tc4_s2
    vle16.v v16, (t1)
    li a0, 0x0001
    SAVE_CSRS
    vnmsac.vx v8, a0, v16
    SET_TEST_NUM 4
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 8
    CHECK_CSRS_UNCHANGED

    /* Test 5: vnmsac.vx SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc5_vd
    vle32.v v8, (t1)
    la t1, tc5_s2
    vle32.v v16, (t1)
    li a0, 0x00000002
    SAVE_CSRS
    vnmsac.vx v8, a0, v16
    SET_TEST_NUM 5
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 16
    CHECK_CSRS_UNCHANGED

    /* Test 6: vnmsac.vx SEW=32 accum */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc6_vd
    vle32.v v8, (t1)
    la t1, tc6_s2
    vle32.v v16, (t1)
    li a0, 0x00000001
    SAVE_CSRS
    vnmsac.vx v8, a0, v16
    SET_TEST_NUM 6
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 16
    CHECK_CSRS_UNCHANGED

    /* Test 7: vnmsac.vx SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc7_vd
    vle64.v v8, (t1)
    la t1, tc7_s2
    vle64.v v16, (t1)
    li a0, 0x0000000000000002
    SAVE_CSRS
    vnmsac.vx v8, a0, v16
    SET_TEST_NUM 7
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 32
    CHECK_CSRS_UNCHANGED

    /* Test 8: vnmsac.vx SEW=64 accum */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc8_vd
    vle64.v v8, (t1)
    la t1, tc8_s2
    vle64.v v16, (t1)
    li a0, 0x0000000000000001
    SAVE_CSRS
    vnmsac.vx v8, a0, v16
    SET_TEST_NUM 8
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_vd:
    .byte 0x00, 0x00, 0x00, 0x00
tc1_s2:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .byte 0xec, 0xd8, 0xc4, 0xb0
.align 1
tc2_vd:
    .byte 0x64, 0xc8, 0x2c, 0x90
tc2_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc2_exp:
    .byte 0x63, 0xc6, 0x29, 0x8c
.align 1
tc3_vd:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc3_s2:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc3_exp:
    .half 0xffec, 0xffd8, 0xffc4, 0xffb0
.align 1
tc4_vd:
    .half 0x0064, 0x00c8, 0x012c, 0x0190
tc4_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc4_exp:
    .half 0x0063, 0x00c6, 0x0129, 0x018c
.align 2
tc5_vd:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc5_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc5_exp:
    .word 0xffffffec, 0xffffffd8, 0xffffffc4, 0xffffffb0
.align 2
tc6_vd:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
tc6_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc6_exp:
    .word 0x00000063, 0x000000c6, 0x00000129, 0x0000018c
.align 3
tc7_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc7_s2:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc7_exp:
    .dword 0xffffffffffffffec, 0xffffffffffffffd8, 0xffffffffffffffc4, 0xffffffffffffffb0
.align 3
tc8_vd:
    .dword 0x0000000000000064, 0x00000000000000c8, 0x000000000000012c, 0x0000000000000190
tc8_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc8_exp:
    .dword 0x0000000000000063, 0x00000000000000c6, 0x0000000000000129, 0x000000000000018c

.align 4
result_buf:  .space 256
witness_buf: .space 256

