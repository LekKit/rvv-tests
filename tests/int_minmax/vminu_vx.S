/* Auto-generated test for vminu.vx
 * Integer vminu (vx form)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vminu.vx e8 basic: wrong result
 *     2 = vminu.vx e8 basic: witness changed
 *     3 = vminu.vx e8 basic: CSR side-effect
 *     4 = vminu.vx e8 zero_sc: wrong result
 *     5 = vminu.vx e8 zero_sc: witness changed
 *     6 = vminu.vx e8 zero_sc: CSR side-effect
 *     7 = vminu.vx e8 max_sc: wrong result
 *     8 = vminu.vx e8 max_sc: witness changed
 *     9 = vminu.vx e8 max_sc: CSR side-effect
 *    10 = vminu.vx e8 max_vec: wrong result
 *    11 = vminu.vx e8 max_vec: witness changed
 *    12 = vminu.vx e8 max_vec: CSR side-effect
 *    13 = vminu.vx e8 signed_sc: wrong result
 *    14 = vminu.vx e8 signed_sc: witness changed
 *    15 = vminu.vx e8 signed_sc: CSR side-effect
 *    16 = vminu.vx e16 basic: wrong result
 *    17 = vminu.vx e16 basic: witness changed
 *    18 = vminu.vx e16 basic: CSR side-effect
 *    19 = vminu.vx e16 zero_sc: wrong result
 *    20 = vminu.vx e16 zero_sc: witness changed
 *    21 = vminu.vx e16 zero_sc: CSR side-effect
 *    22 = vminu.vx e16 max_sc: wrong result
 *    23 = vminu.vx e16 max_sc: witness changed
 *    24 = vminu.vx e16 max_sc: CSR side-effect
 *    25 = vminu.vx e16 max_vec: wrong result
 *    26 = vminu.vx e16 max_vec: witness changed
 *    27 = vminu.vx e16 max_vec: CSR side-effect
 *    28 = vminu.vx e16 signed_sc: wrong result
 *    29 = vminu.vx e16 signed_sc: witness changed
 *    30 = vminu.vx e16 signed_sc: CSR side-effect
 *    31 = vminu.vx e32 basic: wrong result
 *    32 = vminu.vx e32 basic: witness changed
 *    33 = vminu.vx e32 basic: CSR side-effect
 *    34 = vminu.vx e32 zero_sc: wrong result
 *    35 = vminu.vx e32 zero_sc: witness changed
 *    36 = vminu.vx e32 zero_sc: CSR side-effect
 *    37 = vminu.vx e32 max_sc: wrong result
 *    38 = vminu.vx e32 max_sc: witness changed
 *    39 = vminu.vx e32 max_sc: CSR side-effect
 *    40 = vminu.vx e32 max_vec: wrong result
 *    41 = vminu.vx e32 max_vec: witness changed
 *    42 = vminu.vx e32 max_vec: CSR side-effect
 *    43 = vminu.vx e32 signed_sc: wrong result
 *    44 = vminu.vx e32 signed_sc: witness changed
 *    45 = vminu.vx e32 signed_sc: CSR side-effect
 *    46 = vminu.vx e32 basic overlap(vd=vs2): wrong result
 *    47 = vminu.vx e32 basic overlap(vd=vs2): witness changed
 *    48 = vminu.vx e32 basic overlap(vd=vs2): CSR side-effect
 *    49 = vminu.vx e32 basic masked: wrong result
 *    50 = vminu.vx e32 basic masked: witness changed
 *    51 = vminu.vx e32 basic masked: CSR side-effect
 *    52 = vminu.vx e64 basic: wrong result
 *    53 = vminu.vx e64 basic: witness changed
 *    54 = vminu.vx e64 basic: CSR side-effect
 *    55 = vminu.vx e64 zero_sc: wrong result
 *    56 = vminu.vx e64 zero_sc: witness changed
 *    57 = vminu.vx e64 zero_sc: CSR side-effect
 *    58 = vminu.vx e64 max_sc: wrong result
 *    59 = vminu.vx e64 max_sc: witness changed
 *    60 = vminu.vx e64 max_sc: CSR side-effect
 *    61 = vminu.vx e64 max_vec: wrong result
 *    62 = vminu.vx e64 max_vec: witness changed
 *    63 = vminu.vx e64 max_vec: CSR side-effect
 *    64 = vminu.vx e64 signed_sc: wrong result
 *    65 = vminu.vx e64 signed_sc: witness changed
 *    66 = vminu.vx e64 signed_sc: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vminu.vx SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    li a0, 0x0a
    la t1, tc1_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 4-6: vminu.vx SEW=8 zero_sc */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_s2
    vle8.v v16, (t1)
    li a0, 0x00
    la t1, tc4_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 7-9: vminu.vx SEW=8 max_sc */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_s2
    vle8.v v16, (t1)
    li a0, 0xff
    la t1, tc7_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 10-12: vminu.vx SEW=8 max_vec */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc10_s2
    vle8.v v16, (t1)
    li a0, 0x01
    la t1, tc10_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 13-15: vminu.vx SEW=8 signed_sc */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc13_s2
    vle8.v v16, (t1)
    li a0, 0x80
    la t1, tc13_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 16-18: vminu.vx SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc16_s2
    vle16.v v16, (t1)
    li a0, 0x000a
    la t1, tc16_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 19-21: vminu.vx SEW=16 zero_sc */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc19_s2
    vle16.v v16, (t1)
    li a0, 0x0000
    la t1, tc19_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 22-24: vminu.vx SEW=16 max_sc */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc22_s2
    vle16.v v16, (t1)
    li a0, 0xffff
    la t1, tc22_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 25-27: vminu.vx SEW=16 max_vec */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc25_s2
    vle16.v v16, (t1)
    li a0, 0x0001
    la t1, tc25_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 28-30: vminu.vx SEW=16 signed_sc */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc28_s2
    vle16.v v16, (t1)
    li a0, 0x8000
    la t1, tc28_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 31-33: vminu.vx SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc31_s2
    vle32.v v16, (t1)
    li a0, 0x0000000a
    la t1, tc31_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 34-36: vminu.vx SEW=32 zero_sc */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc34_s2
    vle32.v v16, (t1)
    li a0, 0x00000000
    la t1, tc34_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 37-39: vminu.vx SEW=32 max_sc */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc37_s2
    vle32.v v16, (t1)
    li a0, 0xffffffff
    la t1, tc37_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 40-42: vminu.vx SEW=32 max_vec */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc40_s2
    vle32.v v16, (t1)
    li a0, 0x00000001
    la t1, tc40_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 43-45: vminu.vx SEW=32 signed_sc */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc43_s2
    vle32.v v16, (t1)
    li a0, 0x80000000
    la t1, tc43_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 46-48: vminu.vx SEW=32 basic (vd=vs2) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc46_s2
    vle32.v v8, (t1)
    li a0, 0x0000000a
    la t1, tc46_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v8, a0
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

    /* Test 49-51: vminu.vx SEW=32 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc49_s2
    vle32.v v16, (t1)
    li a0, 0x0000000a
    la t1, tc49_vd
    vle32.v v8, (t1)
    la t1, tc49_mask
    vlm.v v0, (t1)
    la t1, tc49_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0, v0.t
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

    /* Test 52-54: vminu.vx SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc52_s2
    vle64.v v16, (t1)
    li a0, 0x000000000000000a
    la t1, tc52_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 55-57: vminu.vx SEW=64 zero_sc */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc55_s2
    vle64.v v16, (t1)
    li a0, 0x0000000000000000
    la t1, tc55_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 58-60: vminu.vx SEW=64 max_sc */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc58_s2
    vle64.v v16, (t1)
    li a0, 0xffffffffffffffff
    la t1, tc58_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 61-63: vminu.vx SEW=64 max_vec */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc61_s2
    vle64.v v16, (t1)
    li a0, 0x0000000000000001
    la t1, tc61_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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

    /* Test 64-66: vminu.vx SEW=64 signed_sc */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc64_s2
    vle64.v v16, (t1)
    li a0, 0x8000000000000000
    la t1, tc64_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vminu.vx v8, v16, a0
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
    .byte 0x00, 0x00, 0x00, 0x00
tc4_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc7_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc7_exp:
    .byte 0x01, 0x02, 0x03, 0x04
tc7_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc10_s2:
    .byte 0xff, 0xff, 0xff, 0xff
tc10_exp:
    .byte 0x01, 0x01, 0x01, 0x01
tc10_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc13_s2:
    .byte 0x00, 0x01, 0x80, 0xff
tc13_exp:
    .byte 0x00, 0x01, 0x80, 0x80
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
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc19_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc22_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc22_exp:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc22_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc25_s2:
    .half 0xffff, 0xffff, 0xffff, 0xffff
tc25_exp:
    .half 0x0001, 0x0001, 0x0001, 0x0001
tc25_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc28_s2:
    .half 0x0000, 0x0001, 0x8000, 0xffff
tc28_exp:
    .half 0x0000, 0x0001, 0x8000, 0x8000
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
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc34_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc37_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc37_exp:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc37_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc40_s2:
    .word 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff
tc40_exp:
    .word 0x00000001, 0x00000001, 0x00000001, 0x00000001
tc40_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc43_s2:
    .word 0x00000000, 0x00000001, 0x80000000, 0xffffffff
tc43_exp:
    .word 0x00000000, 0x00000001, 0x80000000, 0x80000000
tc43_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc46_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc46_exp:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
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
    .word 0x0000dead, 0x00000002, 0x0000dead, 0x00000004
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
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc55_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc58_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc58_exp:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc58_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc61_s2:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
tc61_exp:
    .dword 0x0000000000000001, 0x0000000000000001, 0x0000000000000001, 0x0000000000000001
tc61_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc64_s2:
    .dword 0x0000000000000000, 0x0000000000000001, 0x8000000000000000, 0xffffffffffffffff
tc64_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0x8000000000000000, 0x8000000000000000
tc64_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

