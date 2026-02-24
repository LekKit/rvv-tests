/* Auto-generated test for vredor.vs
 * Reduction vredor.vs
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vredor.vs e8 basic: wrong result
 *     2 = vredor.vs e8 basic: witness changed
 *     3 = vredor.vs e8 basic: CSR side-effect
 *     4 = vredor.vs e8 init: wrong result
 *     5 = vredor.vs e8 init: witness changed
 *     6 = vredor.vs e8 init: CSR side-effect
 *     7 = vredor.vs e8 max: wrong result
 *     8 = vredor.vs e8 max: witness changed
 *     9 = vredor.vs e8 max: CSR side-effect
 *    10 = vredor.vs e8 basic masked: wrong result
 *    11 = vredor.vs e8 basic masked: witness changed
 *    12 = vredor.vs e8 basic masked: CSR side-effect
 *    13 = vredor.vs e16 basic: wrong result
 *    14 = vredor.vs e16 basic: witness changed
 *    15 = vredor.vs e16 basic: CSR side-effect
 *    16 = vredor.vs e16 init: wrong result
 *    17 = vredor.vs e16 init: witness changed
 *    18 = vredor.vs e16 init: CSR side-effect
 *    19 = vredor.vs e16 max: wrong result
 *    20 = vredor.vs e16 max: witness changed
 *    21 = vredor.vs e16 max: CSR side-effect
 *    22 = vredor.vs e16 basic masked: wrong result
 *    23 = vredor.vs e16 basic masked: witness changed
 *    24 = vredor.vs e16 basic masked: CSR side-effect
 *    25 = vredor.vs e32 basic: wrong result
 *    26 = vredor.vs e32 basic: witness changed
 *    27 = vredor.vs e32 basic: CSR side-effect
 *    28 = vredor.vs e32 init: wrong result
 *    29 = vredor.vs e32 init: witness changed
 *    30 = vredor.vs e32 init: CSR side-effect
 *    31 = vredor.vs e32 max: wrong result
 *    32 = vredor.vs e32 max: witness changed
 *    33 = vredor.vs e32 max: CSR side-effect
 *    34 = vredor.vs e32 basic masked: wrong result
 *    35 = vredor.vs e32 basic masked: witness changed
 *    36 = vredor.vs e32 basic masked: CSR side-effect
 *    37 = vredor.vs e64 basic: wrong result
 *    38 = vredor.vs e64 basic: witness changed
 *    39 = vredor.vs e64 basic: CSR side-effect
 *    40 = vredor.vs e64 init: wrong result
 *    41 = vredor.vs e64 init: witness changed
 *    42 = vredor.vs e64 init: CSR side-effect
 *    43 = vredor.vs e64 max: wrong result
 *    44 = vredor.vs e64 max: witness changed
 *    45 = vredor.vs e64 max: CSR side-effect
 *    46 = vredor.vs e64 basic masked: wrong result
 *    47 = vredor.vs e64 basic masked: witness changed
 *    48 = vredor.vs e64 basic masked: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vredor.vs SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_vec
    vle8.v v16, (t1)
    la t1, tc1_init
    vle8.v v20, (t1)
    la t1, tc1_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 1
    SET_TEST_NUM 2
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 4
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED

    /* Test 4-6: vredor.vs SEW=8 init */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_vec
    vle8.v v16, (t1)
    la t1, tc4_init
    vle8.v v20, (t1)
    la t1, tc4_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 1
    SET_TEST_NUM 5
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 4
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-9: vredor.vs SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_vec
    vle8.v v16, (t1)
    la t1, tc7_init
    vle8.v v20, (t1)
    la t1, tc7_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 1
    SET_TEST_NUM 8
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 4
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED

    /* Test 10-12: vredor.vs SEW=8 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc10_vec
    vle8.v v16, (t1)
    la t1, tc10_init
    vle8.v v20, (t1)
    la t1, tc10_mask
    vlm.v v0, (t1)
    la t1, tc10_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20, v0.t
    SET_TEST_NUM 10
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 1
    SET_TEST_NUM 11
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 4
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13-15: vredor.vs SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc13_vec
    vle16.v v16, (t1)
    la t1, tc13_init
    vle16.v v20, (t1)
    la t1, tc13_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 13
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 2
    SET_TEST_NUM 14
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 8
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED

    /* Test 16-18: vredor.vs SEW=16 init */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc16_vec
    vle16.v v16, (t1)
    la t1, tc16_init
    vle16.v v20, (t1)
    la t1, tc16_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 16
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 2
    SET_TEST_NUM 17
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc16_w, 8
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    /* Test 19-21: vredor.vs SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc19_vec
    vle16.v v16, (t1)
    la t1, tc19_init
    vle16.v v20, (t1)
    la t1, tc19_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 19
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 2
    SET_TEST_NUM 20
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc19_w, 8
    SET_TEST_NUM 21
    CHECK_CSRS_UNCHANGED

    /* Test 22-24: vredor.vs SEW=16 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc22_vec
    vle16.v v16, (t1)
    la t1, tc22_init
    vle16.v v20, (t1)
    la t1, tc22_mask
    vlm.v v0, (t1)
    la t1, tc22_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20, v0.t
    SET_TEST_NUM 22
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc22_exp, 2
    SET_TEST_NUM 23
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc22_w, 8
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    /* Test 25-27: vredor.vs SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc25_vec
    vle32.v v16, (t1)
    la t1, tc25_init
    vle32.v v20, (t1)
    la t1, tc25_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 25
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc25_exp, 4
    SET_TEST_NUM 26
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc25_w, 16
    SET_TEST_NUM 27
    CHECK_CSRS_UNCHANGED

    /* Test 28-30: vredor.vs SEW=32 init */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc28_vec
    vle32.v v16, (t1)
    la t1, tc28_init
    vle32.v v20, (t1)
    la t1, tc28_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 28
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc28_exp, 4
    SET_TEST_NUM 29
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc28_w, 16
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED

    /* Test 31-33: vredor.vs SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc31_vec
    vle32.v v16, (t1)
    la t1, tc31_init
    vle32.v v20, (t1)
    la t1, tc31_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 31
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc31_exp, 4
    SET_TEST_NUM 32
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc31_w, 16
    SET_TEST_NUM 33
    CHECK_CSRS_UNCHANGED

    /* Test 34-36: vredor.vs SEW=32 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc34_vec
    vle32.v v16, (t1)
    la t1, tc34_init
    vle32.v v20, (t1)
    la t1, tc34_mask
    vlm.v v0, (t1)
    la t1, tc34_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20, v0.t
    SET_TEST_NUM 34
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc34_exp, 4
    SET_TEST_NUM 35
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc34_w, 16
    SET_TEST_NUM 36
    CHECK_CSRS_UNCHANGED

    /* Test 37-39: vredor.vs SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc37_vec
    vle64.v v16, (t1)
    la t1, tc37_init
    vle64.v v20, (t1)
    la t1, tc37_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 37
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc37_exp, 8
    SET_TEST_NUM 38
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc37_w, 32
    SET_TEST_NUM 39
    CHECK_CSRS_UNCHANGED

    /* Test 40-42: vredor.vs SEW=64 init */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc40_vec
    vle64.v v16, (t1)
    la t1, tc40_init
    vle64.v v20, (t1)
    la t1, tc40_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 40
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc40_exp, 8
    SET_TEST_NUM 41
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc40_w, 32
    SET_TEST_NUM 42
    CHECK_CSRS_UNCHANGED

    /* Test 43-45: vredor.vs SEW=64 max */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc43_vec
    vle64.v v16, (t1)
    la t1, tc43_init
    vle64.v v20, (t1)
    la t1, tc43_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20
    SET_TEST_NUM 43
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc43_exp, 8
    SET_TEST_NUM 44
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc43_w, 32
    SET_TEST_NUM 45
    CHECK_CSRS_UNCHANGED

    /* Test 46-48: vredor.vs SEW=64 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc46_vec
    vle64.v v16, (t1)
    la t1, tc46_init
    vle64.v v20, (t1)
    la t1, tc46_mask
    vlm.v v0, (t1)
    la t1, tc46_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vredor.vs v8, v16, v20, v0.t
    SET_TEST_NUM 46
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc46_exp, 8
    SET_TEST_NUM 47
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc46_w, 32
    SET_TEST_NUM 48
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_vec:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_init:
    .byte 0x00, 0x00, 0x00, 0x00
tc1_exp:
    .byte 0x07, 0x00, 0x00, 0x00
tc1_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc4_vec:
    .byte 0x01, 0x02, 0x03, 0x04
tc4_init:
    .byte 0x64, 0x00, 0x00, 0x00
tc4_exp:
    .byte 0x67, 0x00, 0x00, 0x00
tc4_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc7_vec:
    .byte 0xff, 0xfe, 0xfd, 0x00
tc7_init:
    .byte 0x00, 0x00, 0x00, 0x00
tc7_exp:
    .byte 0xff, 0x00, 0x00, 0x00
tc7_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc10_mask:
    .byte 5
.align 1
tc10_vec:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc10_init:
    .byte 0x64, 0x00, 0x00, 0x00
tc10_exp:
    .byte 0x7e, 0x00, 0x00, 0x00
tc10_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc13_vec:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc13_init:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc13_exp:
    .half 0x0007, 0x0000, 0x0000, 0x0000
tc13_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc16_vec:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc16_init:
    .half 0x0064, 0x0000, 0x0000, 0x0000
tc16_exp:
    .half 0x0067, 0x0000, 0x0000, 0x0000
tc16_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc19_vec:
    .half 0xffff, 0xfffe, 0xfffd, 0x0000
tc19_init:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc19_exp:
    .half 0xffff, 0x0000, 0x0000, 0x0000
tc19_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc22_mask:
    .byte 5
.align 1
tc22_vec:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc22_init:
    .half 0x0064, 0x0000, 0x0000, 0x0000
tc22_exp:
    .half 0x007e, 0x0000, 0x0000, 0x0000
tc22_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 2
tc25_vec:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc25_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc25_exp:
    .word 0x00000007, 0x00000000, 0x00000000, 0x00000000
tc25_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc28_vec:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc28_init:
    .word 0x00000064, 0x00000000, 0x00000000, 0x00000000
tc28_exp:
    .word 0x00000067, 0x00000000, 0x00000000, 0x00000000
tc28_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc31_vec:
    .word 0xffffffff, 0xfffffffe, 0xfffffffd, 0x00000000
tc31_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc31_exp:
    .word 0xffffffff, 0x00000000, 0x00000000, 0x00000000
tc31_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 1
tc34_mask:
    .byte 5
.align 2
tc34_vec:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc34_init:
    .word 0x00000064, 0x00000000, 0x00000000, 0x00000000
tc34_exp:
    .word 0x0000007e, 0x00000000, 0x00000000, 0x00000000
tc34_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc37_vec:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc37_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc37_exp:
    .dword 0x0000000000000007, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc37_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc40_vec:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc40_init:
    .dword 0x0000000000000064, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc40_exp:
    .dword 0x0000000000000067, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc40_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc43_vec:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0xfffffffffffffffd, 0x0000000000000000
tc43_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc43_exp:
    .dword 0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc43_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 1
tc46_mask:
    .byte 5
.align 3
tc46_vec:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc46_init:
    .dword 0x0000000000000064, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc46_exp:
    .dword 0x000000000000007e, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc46_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

