/* Auto-generated test for vrgather.vv
 * Gather elements by index vector
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vrgather.vv e8 identity: result
 *     2 = vrgather.vv e8 reverse: result
 *     3 = vrgather.vv e8 splat: result
 *     4 = vrgather.vv e8 oob: result
 *     5 = vrgather.vv e16 identity: result
 *     6 = vrgather.vv e16 reverse: result
 *     7 = vrgather.vv e16 splat: result
 *     8 = vrgather.vv e16 oob: result
 *     9 = vrgather.vv e32 identity: result
 *    10 = vrgather.vv e32 reverse: result
 *    11 = vrgather.vv e32 splat: result
 *    12 = vrgather.vv e32 oob: result
 *    13 = vrgather.vv e64 identity: result
 *    14 = vrgather.vv e64 reverse: result
 *    15 = vrgather.vv e64 splat: result
 *    16 = vrgather.vv e64 oob: result
 *    17 = vrgather.vv e32 masked: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    vle8.v v16, (t1)
    la t1, tc1_idx
    vle8.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_src
    vle8.v v16, (t1)
    la t1, tc2_idx
    vle8.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_src
    vle8.v v16, (t1)
    la t1, tc3_idx
    vle8.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_src
    vle8.v v16, (t1)
    la t1, tc4_idx
    vle8.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_src
    vle16.v v16, (t1)
    la t1, tc5_idx
    vle16.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 5
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_src
    vle16.v v16, (t1)
    la t1, tc6_idx
    vle16.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 6
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_src
    vle16.v v16, (t1)
    la t1, tc7_idx
    vle16.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc8_src
    vle16.v v16, (t1)
    la t1, tc8_idx
    vle16.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 8
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_src
    vle32.v v16, (t1)
    la t1, tc9_idx
    vle32.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 9
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc10_src
    vle32.v v16, (t1)
    la t1, tc10_idx
    vle32.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 10
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc11_src
    vle32.v v16, (t1)
    la t1, tc11_idx
    vle32.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 11
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc12_src
    vle32.v v16, (t1)
    la t1, tc12_idx
    vle32.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 12
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc13_src
    vle64.v v16, (t1)
    la t1, tc13_idx
    vle64.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 13
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc14_src
    vle64.v v16, (t1)
    la t1, tc14_idx
    vle64.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 14
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc14_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc15_src
    vle64.v v16, (t1)
    la t1, tc15_idx
    vle64.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 15
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc15_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_src
    vle64.v v16, (t1)
    la t1, tc16_idx
    vle64.v v20, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    CHECK_CSRS_UNCHANGED

    /* Masked gather: mask=0b0101, tu/mu */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc17_vd
    vle32.v v8, (t1)
    la t1, tc17_src
    vle32.v v16, (t1)
    la t1, tc17_idx
    vle32.v v20, (t1)
    la t1, tc17_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vrgather.vv v8, v16, v20, v0.t
    SET_TEST_NUM 17
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc17_exp, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_idx:
    .byte 0x00, 0x01, 0x02, 0x03
tc1_exp:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc2_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc2_idx:
    .byte 0x03, 0x02, 0x01, 0x00
tc2_exp:
    .byte 0x28, 0x1e, 0x14, 0x0a
.align 1
tc3_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc3_idx:
    .byte 0x00, 0x00, 0x00, 0x00
tc3_exp:
    .byte 0x0a, 0x0a, 0x0a, 0x0a
.align 1
tc4_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc4_idx:
    .byte 0x00, 0x01, 0x02, 0xff
tc4_exp:
    .byte 0x0a, 0x14, 0x1e, 0x00
.align 1
tc5_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc5_idx:
    .half 0x0000, 0x0001, 0x0002, 0x0003
tc5_exp:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc6_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc6_idx:
    .half 0x0003, 0x0002, 0x0001, 0x0000
tc6_exp:
    .half 0x0028, 0x001e, 0x0014, 0x000a
.align 1
tc7_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc7_idx:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc7_exp:
    .half 0x000a, 0x000a, 0x000a, 0x000a
.align 1
tc8_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc8_idx:
    .half 0x0000, 0x0001, 0x0002, 0xffff
tc8_exp:
    .half 0x000a, 0x0014, 0x001e, 0x0000
.align 2
tc9_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc9_idx:
    .word 0x00000000, 0x00000001, 0x00000002, 0x00000003
tc9_exp:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 2
tc10_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc10_idx:
    .word 0x00000003, 0x00000002, 0x00000001, 0x00000000
tc10_exp:
    .word 0x00000028, 0x0000001e, 0x00000014, 0x0000000a
.align 2
tc11_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc11_idx:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc11_exp:
    .word 0x0000000a, 0x0000000a, 0x0000000a, 0x0000000a
.align 2
tc12_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc12_idx:
    .word 0x00000000, 0x00000001, 0x00000002, 0xffffffff
tc12_exp:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000000
.align 3
tc13_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc13_idx:
    .dword 0x0000000000000000, 0x0000000000000001, 0x0000000000000002, 0x0000000000000003
tc13_exp:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 3
tc14_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc14_idx:
    .dword 0x0000000000000003, 0x0000000000000002, 0x0000000000000001, 0x0000000000000000
tc14_exp:
    .dword 0x0000000000000028, 0x000000000000001e, 0x0000000000000014, 0x000000000000000a
.align 3
tc15_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc15_idx:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc15_exp:
    .dword 0x000000000000000a, 0x000000000000000a, 0x000000000000000a, 0x000000000000000a
.align 3
tc16_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc16_idx:
    .dword 0x0000000000000000, 0x0000000000000001, 0x0000000000000002, 0xffffffffffffffff
tc16_exp:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000000
.align 1
tc17_mask:
    .byte 5
.align 2
tc17_vd:
    .word 0x000000aa, 0x000000bb, 0x000000cc, 0x000000dd
tc17_src:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
tc17_idx:
    .word 0x00000003, 0x00000002, 0x00000001, 0x00000000
tc17_exp:
    .word 0x00000190, 0x000000bb, 0x000000c8, 0x000000dd

.align 4
result_buf:  .space 256
witness_buf: .space 256

