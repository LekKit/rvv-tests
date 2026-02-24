/* Auto-generated test for vrgatherei16.vv
 * Gather with 16-bit index vector
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vrgatherei16 e8 identity: result
 *     2 = vrgatherei16 e8 reverse: result
 *     3 = vrgatherei16 e8 splat: result
 *     4 = vrgatherei16 e16 identity: result
 *     5 = vrgatherei16 e16 reverse: result
 *     6 = vrgatherei16 e16 splat: result
 *     7 = vrgatherei16 e32 identity: result
 *     8 = vrgatherei16 e32 reverse: result
 *     9 = vrgatherei16 e32 splat: result
 *    10 = vrgatherei16 e64 identity: result
 *    11 = vrgatherei16 e64 reverse: result
 *    12 = vrgatherei16 e64 splat: result
 *    13 = vrgatherei16 e32 masked: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc2_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 2
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc3_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc4_src
    vle16.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc4_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_src
    vle16.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc5_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 5
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_src
    vle16.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc6_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 6
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc6_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc7_src
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc7_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc8_src
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc8_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 8
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc8_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc9_src
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc9_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 9
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc9_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc10_src
    vle64.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc10_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 10
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc11_src
    vle64.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc11_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 11
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc11_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc12_src
    vle64.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc12_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e64, m1, tu, mu
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20
    SET_TEST_NUM 12
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc12_exp, 32
    CHECK_CSRS_UNCHANGED

    /* Masked vrgatherei16.vv */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_vd
    vle32.v v8, (t1)
    la t1, tc13_src
    vle32.v v16, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc13_idx
    vle16.v v20, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vrgatherei16.vv v8, v16, v20, v0.t
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc1_idx:
    .half 0x0000, 0x0001, 0x0002, 0x0003
.align 1
tc1_exp:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc2_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc2_idx:
    .half 0x0003, 0x0002, 0x0001, 0x0000
.align 1
tc2_exp:
    .byte 0x28, 0x1e, 0x14, 0x0a
.align 1
tc3_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc3_idx:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 1
tc3_exp:
    .byte 0x0a, 0x0a, 0x0a, 0x0a
.align 1
tc4_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc4_idx:
    .half 0x0000, 0x0001, 0x0002, 0x0003
.align 1
tc4_exp:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc5_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc5_idx:
    .half 0x0003, 0x0002, 0x0001, 0x0000
.align 1
tc5_exp:
    .half 0x0028, 0x001e, 0x0014, 0x000a
.align 1
tc6_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 1
tc6_idx:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 1
tc6_exp:
    .half 0x000a, 0x000a, 0x000a, 0x000a
.align 2
tc7_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 1
tc7_idx:
    .half 0x0000, 0x0001, 0x0002, 0x0003
.align 2
tc7_exp:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 2
tc8_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 1
tc8_idx:
    .half 0x0003, 0x0002, 0x0001, 0x0000
.align 2
tc8_exp:
    .word 0x00000028, 0x0000001e, 0x00000014, 0x0000000a
.align 2
tc9_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 1
tc9_idx:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 2
tc9_exp:
    .word 0x0000000a, 0x0000000a, 0x0000000a, 0x0000000a
.align 3
tc10_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 1
tc10_idx:
    .half 0x0000, 0x0001, 0x0002, 0x0003
.align 3
tc10_exp:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 3
tc11_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 1
tc11_idx:
    .half 0x0003, 0x0002, 0x0001, 0x0000
.align 3
tc11_exp:
    .dword 0x0000000000000028, 0x000000000000001e, 0x0000000000000014, 0x000000000000000a
.align 3
tc12_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
.align 1
tc12_idx:
    .half 0x0000, 0x0000, 0x0000, 0x0000
.align 3
tc12_exp:
    .dword 0x000000000000000a, 0x000000000000000a, 0x000000000000000a, 0x000000000000000a
.align 1
tc13_mask:
    .byte 5
.align 2
tc13_vd:
    .word 0x000000aa, 0x000000bb, 0x000000cc, 0x000000dd
tc13_src:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
.align 1
tc13_idx:
    .half 0x0003, 0x0002, 0x0001, 0x0000
.align 2
tc13_exp:
    .word 0x00000190, 0x000000bb, 0x000000c8, 0x000000dd

.align 4
result_buf:  .space 256
witness_buf: .space 256

