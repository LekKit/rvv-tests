/* Auto-generated test for vslide1up.vx
 * Slide1 up
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vslide1up.vx e8: result
 *     2 = vslide1up.vx e16: result
 *     3 = vslide1up.vx e32: result
 *     4 = vslide1up.vx e64: result
 *     5 = vslide1up.vx e32 masked: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    vle8.v v16, (t1)
    li a0, 99
    SAVE_CSRS
    vslide1up.vx v8, v16, a0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_src
    vle16.v v16, (t1)
    li a0, 99
    SAVE_CSRS
    vslide1up.vx v8, v16, a0
    SET_TEST_NUM 2
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_src
    vle32.v v16, (t1)
    li a0, 99
    SAVE_CSRS
    vslide1up.vx v8, v16, a0
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc4_src
    vle64.v v16, (t1)
    li a0, 99
    SAVE_CSRS
    vslide1up.vx v8, v16, a0
    SET_TEST_NUM 4
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 32
    CHECK_CSRS_UNCHANGED

    /* Masked slide1 up: mask=0b0101 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc5_src
    vle32.v v16, (t1)
    la t1, tc5_vd
    vle32.v v8, (t1)
    la t1, tc5_mask
    vlm.v v0, (t1)
    li a0, 99
    SAVE_CSRS
    vslide1up.vx v8, v16, a0, v0.t
    SET_TEST_NUM 5
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .byte 0x63, 0x0a, 0x14, 0x1e
.align 1
tc2_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc2_exp:
    .half 0x0063, 0x000a, 0x0014, 0x001e
.align 2
tc3_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc3_exp:
    .word 0x00000063, 0x0000000a, 0x00000014, 0x0000001e
.align 3
tc4_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc4_exp:
    .dword 0x0000000000000063, 0x000000000000000a, 0x0000000000000014, 0x000000000000001e
.align 1
tc5_mask:
    .byte 5
.align 2
tc5_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc5_vd:
    .word 0x000000aa, 0x000000bb, 0x000000cc, 0x000000dd
tc5_exp:
    .word 0x00000063, 0x000000bb, 0x00000014, 0x000000dd

.align 4
result_buf:  .space 256
witness_buf: .space 256

