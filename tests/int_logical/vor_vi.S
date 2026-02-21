/* Auto-generated test for vor.vi
 * Integer vor (vi form)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vor.vi e8 imm0: wrong result
 *     2 = vor.vi e8 imm0: witness changed
 *     3 = vor.vi e8 imm0: CSR side-effect
 *     4 = vor.vi e8 imm1: wrong result
 *     5 = vor.vi e8 imm1: witness changed
 *     6 = vor.vi e8 imm1: CSR side-effect
 *     7 = vor.vi e8 imm_n1: wrong result
 *     8 = vor.vi e8 imm_n1: witness changed
 *     9 = vor.vi e8 imm_n1: CSR side-effect
 *    10 = vor.vi e8 imm15: wrong result
 *    11 = vor.vi e8 imm15: witness changed
 *    12 = vor.vi e8 imm15: CSR side-effect
 *    13 = vor.vi e8 imm_n16: wrong result
 *    14 = vor.vi e8 imm_n16: witness changed
 *    15 = vor.vi e8 imm_n16: CSR side-effect
 *    16 = vor.vi e16 imm0: wrong result
 *    17 = vor.vi e16 imm0: witness changed
 *    18 = vor.vi e16 imm0: CSR side-effect
 *    19 = vor.vi e16 imm1: wrong result
 *    20 = vor.vi e16 imm1: witness changed
 *    21 = vor.vi e16 imm1: CSR side-effect
 *    22 = vor.vi e16 imm_n1: wrong result
 *    23 = vor.vi e16 imm_n1: witness changed
 *    24 = vor.vi e16 imm_n1: CSR side-effect
 *    25 = vor.vi e16 imm15: wrong result
 *    26 = vor.vi e16 imm15: witness changed
 *    27 = vor.vi e16 imm15: CSR side-effect
 *    28 = vor.vi e16 imm_n16: wrong result
 *    29 = vor.vi e16 imm_n16: witness changed
 *    30 = vor.vi e16 imm_n16: CSR side-effect
 *    31 = vor.vi e32 imm0: wrong result
 *    32 = vor.vi e32 imm0: witness changed
 *    33 = vor.vi e32 imm0: CSR side-effect
 *    34 = vor.vi e32 imm1: wrong result
 *    35 = vor.vi e32 imm1: witness changed
 *    36 = vor.vi e32 imm1: CSR side-effect
 *    37 = vor.vi e32 imm_n1: wrong result
 *    38 = vor.vi e32 imm_n1: witness changed
 *    39 = vor.vi e32 imm_n1: CSR side-effect
 *    40 = vor.vi e32 imm15: wrong result
 *    41 = vor.vi e32 imm15: witness changed
 *    42 = vor.vi e32 imm15: CSR side-effect
 *    43 = vor.vi e32 imm_n16: wrong result
 *    44 = vor.vi e32 imm_n16: witness changed
 *    45 = vor.vi e32 imm_n16: CSR side-effect
 *    46 = vor.vi e32 basic overlap(vd=vs2): wrong result
 *    47 = vor.vi e32 basic overlap(vd=vs2): witness changed
 *    48 = vor.vi e32 basic overlap(vd=vs2): CSR side-effect
 *    49 = vor.vi e32 basic masked: wrong result
 *    50 = vor.vi e32 basic masked: witness changed
 *    51 = vor.vi e32 basic masked: CSR side-effect
 *    52 = vor.vi e64 imm0: wrong result
 *    53 = vor.vi e64 imm0: witness changed
 *    54 = vor.vi e64 imm0: CSR side-effect
 *    55 = vor.vi e64 imm1: wrong result
 *    56 = vor.vi e64 imm1: witness changed
 *    57 = vor.vi e64 imm1: CSR side-effect
 *    58 = vor.vi e64 imm_n1: wrong result
 *    59 = vor.vi e64 imm_n1: witness changed
 *    60 = vor.vi e64 imm_n1: CSR side-effect
 *    61 = vor.vi e64 imm15: wrong result
 *    62 = vor.vi e64 imm15: witness changed
 *    63 = vor.vi e64 imm15: CSR side-effect
 *    64 = vor.vi e64 imm_n16: wrong result
 *    65 = vor.vi e64 imm_n16: witness changed
 *    66 = vor.vi e64 imm_n16: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vor.vi SEW=8 imm0 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    SET_TEST_NUM 2
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 4
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED

    /* Test 4-6: vor.vi SEW=8 imm1 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    la t1, tc4_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 1
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 4
    SET_TEST_NUM 5
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 4
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-9: vor.vi SEW=8 imm_n1 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_s2
    vle8.v v16, (t1)
    la t1, tc7_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, -1
    SET_TEST_NUM 7
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 4
    SET_TEST_NUM 8
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 4
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED

    /* Test 10-12: vor.vi SEW=8 imm15 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc10_s2
    vle8.v v16, (t1)
    la t1, tc10_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 15
    SET_TEST_NUM 10
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 4
    SET_TEST_NUM 11
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 4
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13-15: vor.vi SEW=8 imm_n16 */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc13_s2
    vle8.v v16, (t1)
    la t1, tc13_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, -16
    SET_TEST_NUM 13
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 4
    SET_TEST_NUM 14
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 4
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED

    /* Test 16-18: vor.vi SEW=16 imm0 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc16_s2
    vle16.v v16, (t1)
    la t1, tc16_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 0
    SET_TEST_NUM 16
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 8
    SET_TEST_NUM 17
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc16_w, 8
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    /* Test 19-21: vor.vi SEW=16 imm1 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc19_s2
    vle16.v v16, (t1)
    la t1, tc19_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 1
    SET_TEST_NUM 19
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 8
    SET_TEST_NUM 20
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc19_w, 8
    SET_TEST_NUM 21
    CHECK_CSRS_UNCHANGED

    /* Test 22-24: vor.vi SEW=16 imm_n1 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc22_s2
    vle16.v v16, (t1)
    la t1, tc22_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, -1
    SET_TEST_NUM 22
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc22_exp, 8
    SET_TEST_NUM 23
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc22_w, 8
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    /* Test 25-27: vor.vi SEW=16 imm15 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc25_s2
    vle16.v v16, (t1)
    la t1, tc25_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 15
    SET_TEST_NUM 25
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc25_exp, 8
    SET_TEST_NUM 26
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc25_w, 8
    SET_TEST_NUM 27
    CHECK_CSRS_UNCHANGED

    /* Test 28-30: vor.vi SEW=16 imm_n16 */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc28_s2
    vle16.v v16, (t1)
    la t1, tc28_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, -16
    SET_TEST_NUM 28
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc28_exp, 8
    SET_TEST_NUM 29
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc28_w, 8
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED

    /* Test 31-33: vor.vi SEW=32 imm0 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc31_s2
    vle32.v v16, (t1)
    la t1, tc31_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 0
    SET_TEST_NUM 31
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc31_exp, 16
    SET_TEST_NUM 32
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc31_w, 16
    SET_TEST_NUM 33
    CHECK_CSRS_UNCHANGED

    /* Test 34-36: vor.vi SEW=32 imm1 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc34_s2
    vle32.v v16, (t1)
    la t1, tc34_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 1
    SET_TEST_NUM 34
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc34_exp, 16
    SET_TEST_NUM 35
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc34_w, 16
    SET_TEST_NUM 36
    CHECK_CSRS_UNCHANGED

    /* Test 37-39: vor.vi SEW=32 imm_n1 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc37_s2
    vle32.v v16, (t1)
    la t1, tc37_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, -1
    SET_TEST_NUM 37
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc37_exp, 16
    SET_TEST_NUM 38
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc37_w, 16
    SET_TEST_NUM 39
    CHECK_CSRS_UNCHANGED

    /* Test 40-42: vor.vi SEW=32 imm15 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc40_s2
    vle32.v v16, (t1)
    la t1, tc40_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 15
    SET_TEST_NUM 40
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc40_exp, 16
    SET_TEST_NUM 41
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc40_w, 16
    SET_TEST_NUM 42
    CHECK_CSRS_UNCHANGED

    /* Test 43-45: vor.vi SEW=32 imm_n16 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc43_s2
    vle32.v v16, (t1)
    la t1, tc43_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, -16
    SET_TEST_NUM 43
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc43_exp, 16
    SET_TEST_NUM 44
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc43_w, 16
    SET_TEST_NUM 45
    CHECK_CSRS_UNCHANGED

    /* Test 46-48: vor.vi SEW=32 basic (vd=vs2) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc46_s2
    vle32.v v8, (t1)
    la t1, tc46_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v8, 1
    SET_TEST_NUM 46
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc46_exp, 16
    SET_TEST_NUM 47
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc46_w, 16
    SET_TEST_NUM 48
    CHECK_CSRS_UNCHANGED

    /* Test 49-51: vor.vi SEW=32 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc49_s2
    vle32.v v16, (t1)
    la t1, tc49_vd
    vle32.v v8, (t1)
    la t1, tc49_mask
    vlm.v v0, (t1)
    la t1, tc49_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 1, v0.t
    SET_TEST_NUM 49
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc49_exp, 16
    SET_TEST_NUM 50
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc49_w, 16
    SET_TEST_NUM 51
    CHECK_CSRS_UNCHANGED

    /* Test 52-54: vor.vi SEW=64 imm0 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc52_s2
    vle64.v v16, (t1)
    la t1, tc52_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 0
    SET_TEST_NUM 52
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc52_exp, 32
    SET_TEST_NUM 53
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc52_w, 32
    SET_TEST_NUM 54
    CHECK_CSRS_UNCHANGED

    /* Test 55-57: vor.vi SEW=64 imm1 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc55_s2
    vle64.v v16, (t1)
    la t1, tc55_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 1
    SET_TEST_NUM 55
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc55_exp, 32
    SET_TEST_NUM 56
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc55_w, 32
    SET_TEST_NUM 57
    CHECK_CSRS_UNCHANGED

    /* Test 58-60: vor.vi SEW=64 imm_n1 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc58_s2
    vle64.v v16, (t1)
    la t1, tc58_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, -1
    SET_TEST_NUM 58
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc58_exp, 32
    SET_TEST_NUM 59
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc58_w, 32
    SET_TEST_NUM 60
    CHECK_CSRS_UNCHANGED

    /* Test 61-63: vor.vi SEW=64 imm15 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc61_s2
    vle64.v v16, (t1)
    la t1, tc61_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, 15
    SET_TEST_NUM 61
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc61_exp, 32
    SET_TEST_NUM 62
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc61_w, 32
    SET_TEST_NUM 63
    CHECK_CSRS_UNCHANGED

    /* Test 64-66: vor.vi SEW=64 imm_n16 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc64_s2
    vle64.v v16, (t1)
    la t1, tc64_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vor.vi v8, v16, -16
    SET_TEST_NUM 64
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc64_exp, 32
    SET_TEST_NUM 65
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc64_w, 32
    SET_TEST_NUM 66
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_exp:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc4_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc4_exp:
    .byte 0x01, 0x03, 0x03, 0x05
tc4_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc7_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc7_exp:
    .byte 0xff, 0xff, 0xff, 0xff
tc7_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc10_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc10_exp:
    .byte 0x0f, 0x0f, 0x0f, 0x0f
tc10_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc13_s2:
    .byte 0x00, 0xff, 0x01, 0xfe
tc13_exp:
    .byte 0xf0, 0xff, 0xf1, 0xfe
tc13_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc16_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc16_exp:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc16_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc19_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc19_exp:
    .half 0x0001, 0x0003, 0x0003, 0x0005
tc19_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc22_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc22_exp:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc22_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc25_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc25_exp:
    .half 0x000f, 0x000f, 0x000f, 0x000f
tc25_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc28_s2:
    .half 0x0000, 0xffff, 0x0001, 0xfffe
tc28_exp:
    .half 0xfff0, 0xffff, 0xfff1, 0xfffe
tc28_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 2
tc31_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc31_exp:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc31_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc34_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc34_exp:
    .word 0x00000001, 0x00000003, 0x00000003, 0x00000005
tc34_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc37_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc37_exp:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc37_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc40_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc40_exp:
    .word 0x0000000f, 0x0000000f, 0x0000000f, 0x0000000f
tc40_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc43_s2:
    .word 0x00000000, 0xffffffff, 0x00000001, 0xfffffffe
tc43_exp:
    .word 0xfffffff0, 0xffffffff, 0xfffffff1, 0xfffffffe
tc43_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc46_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc46_exp:
    .word 0x00000001, 0x00000003, 0x00000003, 0x00000005
tc46_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 1
tc49_mask:
    .byte 10
.align 2
tc49_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc49_vd:
    .word 0x0000dead, 0x0000dead, 0x0000dead, 0x0000dead
tc49_exp:
    .word 0x0000dead, 0x00000003, 0x0000dead, 0x00000005
tc49_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc52_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc52_exp:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc52_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc55_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc55_exp:
    .dword 0x0000000000000001, 0x0000000000000003, 0x0000000000000003, 0x0000000000000005
tc55_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc58_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc58_exp:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
tc58_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc61_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc61_exp:
    .dword 0x000000000000000f, 0x000000000000000f, 0x000000000000000f, 0x000000000000000f
tc61_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc64_s2:
    .dword 0x0000000000000000, 0xffffffffffffffff, 0x0000000000000001, 0xfffffffffffffffe
tc64_exp:
    .dword 0xfffffffffffffff0, 0xffffffffffffffff, 0xfffffffffffffff1, 0xfffffffffffffffe
tc64_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

