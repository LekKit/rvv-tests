/* Auto-generated test for vcompress.vm
 * Compress active elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vcompress e8 all: result
 *     2 = vcompress e8 none: result
 *     3 = vcompress e8 even: result
 *     4 = vcompress e8 odd: result
 *     5 = vcompress e8 first: result
 *     6 = vcompress e16 all: result
 *     7 = vcompress e16 none: result
 *     8 = vcompress e16 even: result
 *     9 = vcompress e16 odd: result
 *    10 = vcompress e16 first: result
 *    11 = vcompress e32 all: result
 *    12 = vcompress e32 none: result
 *    13 = vcompress e32 even: result
 *    14 = vcompress e32 odd: result
 *    15 = vcompress e32 first: result
 *    16 = vcompress e64 all: result
 *    17 = vcompress e64 none: result
 *    18 = vcompress e64 even: result
 *    19 = vcompress e64 odd: result
 *    20 = vcompress e64 first: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    vle8.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc1_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
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
    la t1, tc2_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
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
    la t1, tc3_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
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
    la t1, tc4_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_src
    vle8.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc5_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 5
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_src
    vle16.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc6_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
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
    la t1, tc7_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
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
    la t1, tc8_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 8
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_src
    vle16.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc9_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 9
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc10_src
    vle16.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc10_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 10
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_src
    vle32.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc11_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
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
    la t1, tc12_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 12
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_src
    vle32.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc13_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc14_src
    vle32.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc14_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 14
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc14_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc15_src
    vle32.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc15_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 15
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_src
    vle64.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc16_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc17_src
    vle64.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc17_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 17
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc18_src
    vle64.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc18_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 18
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc18_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc19_src
    vle64.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc19_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 19
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc20_src
    vle64.v v16, (t1)
    vmv.v.i v8, 0
    la t1, tc20_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vcompress.vm v8, v16, v0
    SET_TEST_NUM 20
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc20_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_mask:
    .byte 15
.align 1
tc1_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc2_mask:
    .byte 0
.align 1
tc2_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc2_exp:
    .byte 0x00, 0x00, 0x00, 0x00
.align 1
tc3_mask:
    .byte 5
.align 1
tc3_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc3_exp:
    .byte 0x0a, 0x1e, 0x00, 0x00
.align 1
tc4_mask:
    .byte 10
.align 1
tc4_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc4_exp:
    .byte 0x14, 0x28, 0x00, 0x00
.align 1
tc5_mask:
    .byte 1
.align 1
tc5_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc5_exp:
    .byte 0x0a, 0x00, 0x00, 0x00
.align 1
tc6_mask:
    .byte 15
.align 1
tc6_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc6_exp:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc7_mask:
    .byte 0
.align 1
tc7_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc7_exp:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 1
tc8_mask:
    .byte 5
.align 1
tc8_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc8_exp:
    .half 0x000a, 0x001e, 0x0000, 0x0000
.align 1
tc9_mask:
    .byte 10
.align 1
tc9_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc9_exp:
    .half 0x0014, 0x0028, 0x0000, 0x0000
.align 1
tc10_mask:
    .byte 1
.align 1
tc10_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc10_exp:
    .half 0x000a, 0x0000, 0x0000, 0x0000
.align 1
tc11_mask:
    .byte 15
.align 2
tc11_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc11_exp:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 1
tc12_mask:
    .byte 0
.align 2
tc12_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc12_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
.align 1
tc13_mask:
    .byte 5
.align 2
tc13_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc13_exp:
    .word 0x0000000a, 0x0000001e, 0x00000000, 0x00000000
.align 1
tc14_mask:
    .byte 10
.align 2
tc14_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc14_exp:
    .word 0x00000014, 0x00000028, 0x00000000, 0x00000000
.align 1
tc15_mask:
    .byte 1
.align 2
tc15_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc15_exp:
    .word 0x0000000a, 0x00000000, 0x00000000, 0x00000000
.align 1
tc16_mask:
    .byte 15
.align 3
tc16_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc16_exp:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 1
tc17_mask:
    .byte 0
.align 3
tc17_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc17_exp:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
.align 1
tc18_mask:
    .byte 5
.align 3
tc18_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc18_exp:
    .dword 0x000000000000000a, 0x000000000000001e, 0x0000000000000000, 0x0000000000000000
.align 1
tc19_mask:
    .byte 10
.align 3
tc19_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc19_exp:
    .dword 0x0000000000000014, 0x0000000000000028, 0x0000000000000000, 0x0000000000000000
.align 1
tc20_mask:
    .byte 1
.align 3
tc20_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc20_exp:
    .dword 0x000000000000000a, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

