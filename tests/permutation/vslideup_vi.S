/* Auto-generated test for vslideup.vi
 * Slide up by immediate
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vslideup.vi e8 off=0: result
 *     2 = vslideup.vi e8 off=1: result
 *     3 = vslideup.vi e8 off=2: result
 *     4 = vslideup.vi e8 off=3: result
 *     5 = vslideup.vi e16 off=0: result
 *     6 = vslideup.vi e16 off=1: result
 *     7 = vslideup.vi e16 off=2: result
 *     8 = vslideup.vi e16 off=3: result
 *     9 = vslideup.vi e32 off=0: result
 *    10 = vslideup.vi e32 off=1: result
 *    11 = vslideup.vi e32 off=2: result
 *    12 = vslideup.vi e32 off=3: result
 *    13 = vslideup.vi e64 off=0: result
 *    14 = vslideup.vi e64 off=1: result
 *    15 = vslideup.vi e64 off=2: result
 *    16 = vslideup.vi e64 off=3: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    vle8.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_src
    vle8.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 1
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_src
    vle8.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 2
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_src
    vle8.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 3
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_src
    vle16.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 0
    SET_TEST_NUM 5
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_src
    vle16.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 1
    SET_TEST_NUM 6
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_src
    vle16.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 2
    SET_TEST_NUM 7
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc8_src
    vle16.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 3
    SET_TEST_NUM 8
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_src
    vle32.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 0
    SET_TEST_NUM 9
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc10_src
    vle32.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 1
    SET_TEST_NUM 10
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_src
    vle32.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 2
    SET_TEST_NUM 11
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc12_src
    vle32.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 3
    SET_TEST_NUM 12
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc13_src
    vle64.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 0
    SET_TEST_NUM 13
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc14_src
    vle64.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 1
    SET_TEST_NUM 14
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc14_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc15_src
    vle64.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 2
    SET_TEST_NUM 15
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_src
    vle64.v v16, (t1)
    vmv.v.i v8, 0
    SAVE_CSRS
    vslideup.vi v8, v16, 3
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc2_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc2_exp:
    .byte 0x00, 0x0a, 0x14, 0x1e
.align 1
tc3_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc3_exp:
    .byte 0x00, 0x00, 0x0a, 0x14
.align 1
tc4_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc4_exp:
    .byte 0x00, 0x00, 0x00, 0x0a
.align 1
tc5_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc5_exp:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc6_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc6_exp:
    .half 0x0000, 0x000a, 0x0014, 0x001e
.align 1
tc7_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc7_exp:
    .half 0x0000, 0x0000, 0x000a, 0x0014
.align 1
tc8_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc8_exp:
    .half 0x0000, 0x0000, 0x0000, 0x000a
.align 2
tc9_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc9_exp:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 2
tc10_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc10_exp:
    .word 0x00000000, 0x0000000a, 0x00000014, 0x0000001e
.align 2
tc11_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc11_exp:
    .word 0x00000000, 0x00000000, 0x0000000a, 0x00000014
.align 2
tc12_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc12_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x0000000a
.align 3
tc13_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc13_exp:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 3
tc14_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc14_exp:
    .dword 0x0000000000000000, 0x000000000000000a, 0x0000000000000014, 0x000000000000001e
.align 3
tc15_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc15_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x000000000000000a, 0x0000000000000014
.align 3
tc16_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc16_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x000000000000000a

.align 4
result_buf:  .space 256
witness_buf: .space 256

