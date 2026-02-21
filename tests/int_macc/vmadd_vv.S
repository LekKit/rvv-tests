/* Auto-generated test for vmadd.vv
 * Multiply-accumulate vmadd.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmadd.vv e8 basic: wrong result
 *     2 = vmadd.vv e8 basic: witness changed
 *     3 = vmadd.vv e8 basic: CSR side-effect
 *     4 = vmadd.vv e8 accum: wrong result
 *     5 = vmadd.vv e8 accum: witness changed
 *     6 = vmadd.vv e8 accum: CSR side-effect
 *     7 = vmadd.vv e8 max: wrong result
 *     8 = vmadd.vv e8 max: witness changed
 *     9 = vmadd.vv e8 max: CSR side-effect
 *    10 = vmadd.vv e16 basic: wrong result
 *    11 = vmadd.vv e16 basic: witness changed
 *    12 = vmadd.vv e16 basic: CSR side-effect
 *    13 = vmadd.vv e16 accum: wrong result
 *    14 = vmadd.vv e16 accum: witness changed
 *    15 = vmadd.vv e16 accum: CSR side-effect
 *    16 = vmadd.vv e16 max: wrong result
 *    17 = vmadd.vv e16 max: witness changed
 *    18 = vmadd.vv e16 max: CSR side-effect
 *    19 = vmadd.vv e32 basic: wrong result
 *    20 = vmadd.vv e32 basic: witness changed
 *    21 = vmadd.vv e32 basic: CSR side-effect
 *    22 = vmadd.vv e32 accum: wrong result
 *    23 = vmadd.vv e32 accum: witness changed
 *    24 = vmadd.vv e32 accum: CSR side-effect
 *    25 = vmadd.vv e32 max: wrong result
 *    26 = vmadd.vv e32 max: witness changed
 *    27 = vmadd.vv e32 max: CSR side-effect
 *    28 = vmadd.vv e32 basic overlap(vd=vs2): wrong result
 *    29 = vmadd.vv e32 basic overlap(vd=vs2): witness changed
 *    30 = vmadd.vv e32 basic overlap(vd=vs2): CSR side-effect
 *    31 = vmadd.vv e32 basic masked: wrong result
 *    32 = vmadd.vv e32 basic masked: witness changed
 *    33 = vmadd.vv e32 basic masked: CSR side-effect
 *    34 = vmadd.vv e64 basic: wrong result
 *    35 = vmadd.vv e64 basic: witness changed
 *    36 = vmadd.vv e64 basic: CSR side-effect
 *    37 = vmadd.vv e64 accum: wrong result
 *    38 = vmadd.vv e64 accum: witness changed
 *    39 = vmadd.vv e64 accum: CSR side-effect
 *    40 = vmadd.vv e64 max: wrong result
 *    41 = vmadd.vv e64 max: witness changed
 *    42 = vmadd.vv e64 max: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vmadd.vv SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_vd
    vle8.v v8, (t1)
    la t1, tc1_s1
    vle8.v v20, (t1)
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
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

    /* Test 4-6: vmadd.vv SEW=8 accum */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_vd
    vle8.v v8, (t1)
    la t1, tc4_s1
    vle8.v v20, (t1)
    la t1, tc4_s2
    vle8.v v16, (t1)
    la t1, tc4_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
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

    /* Test 7-9: vmadd.vv SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_vd
    vle8.v v8, (t1)
    la t1, tc7_s1
    vle8.v v20, (t1)
    la t1, tc7_s2
    vle8.v v16, (t1)
    la t1, tc7_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
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

    /* Test 10-12: vmadd.vv SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc10_vd
    vle16.v v8, (t1)
    la t1, tc10_s1
    vle16.v v20, (t1)
    la t1, tc10_s2
    vle16.v v16, (t1)
    la t1, tc10_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
    SET_TEST_NUM 10
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 8
    SET_TEST_NUM 11
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 8
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13-15: vmadd.vv SEW=16 accum */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc13_vd
    vle16.v v8, (t1)
    la t1, tc13_s1
    vle16.v v20, (t1)
    la t1, tc13_s2
    vle16.v v16, (t1)
    la t1, tc13_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
    SET_TEST_NUM 13
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 8
    SET_TEST_NUM 14
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 8
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED

    /* Test 16-18: vmadd.vv SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc16_vd
    vle16.v v8, (t1)
    la t1, tc16_s1
    vle16.v v20, (t1)
    la t1, tc16_s2
    vle16.v v16, (t1)
    la t1, tc16_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
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

    /* Test 19-21: vmadd.vv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc19_vd
    vle32.v v8, (t1)
    la t1, tc19_s1
    vle32.v v20, (t1)
    la t1, tc19_s2
    vle32.v v16, (t1)
    la t1, tc19_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
    SET_TEST_NUM 19
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 16
    SET_TEST_NUM 20
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc19_w, 16
    SET_TEST_NUM 21
    CHECK_CSRS_UNCHANGED

    /* Test 22-24: vmadd.vv SEW=32 accum */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc22_vd
    vle32.v v8, (t1)
    la t1, tc22_s1
    vle32.v v20, (t1)
    la t1, tc22_s2
    vle32.v v16, (t1)
    la t1, tc22_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
    SET_TEST_NUM 22
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc22_exp, 16
    SET_TEST_NUM 23
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc22_w, 16
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    /* Test 25-27: vmadd.vv SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc25_vd
    vle32.v v8, (t1)
    la t1, tc25_s1
    vle32.v v20, (t1)
    la t1, tc25_s2
    vle32.v v16, (t1)
    la t1, tc25_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
    SET_TEST_NUM 25
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc25_exp, 16
    SET_TEST_NUM 26
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc25_w, 16
    SET_TEST_NUM 27
    CHECK_CSRS_UNCHANGED

    /* Test 28-30: vmadd.vv SEW=32 basic (vd=vs2) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc28_vd
    vle32.v v8, (t1)
    la t1, tc28_s1
    vle32.v v20, (t1)
    la t1, tc28_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v8
    SET_TEST_NUM 28
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc28_exp, 16
    SET_TEST_NUM 29
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc28_w, 16
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED

    /* Test 31-33: vmadd.vv SEW=32 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc31_vd
    vle32.v v8, (t1)
    la t1, tc31_s1
    vle32.v v20, (t1)
    la t1, tc31_s2
    vle32.v v16, (t1)
    la t1, tc31_mask
    vlm.v v0, (t1)
    la t1, tc31_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16, v0.t
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

    /* Test 34-36: vmadd.vv SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc34_vd
    vle64.v v8, (t1)
    la t1, tc34_s1
    vle64.v v20, (t1)
    la t1, tc34_s2
    vle64.v v16, (t1)
    la t1, tc34_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
    SET_TEST_NUM 34
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc34_exp, 32
    SET_TEST_NUM 35
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc34_w, 32
    SET_TEST_NUM 36
    CHECK_CSRS_UNCHANGED

    /* Test 37-39: vmadd.vv SEW=64 accum */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc37_vd
    vle64.v v8, (t1)
    la t1, tc37_s1
    vle64.v v20, (t1)
    la t1, tc37_s2
    vle64.v v16, (t1)
    la t1, tc37_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
    SET_TEST_NUM 37
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc37_exp, 32
    SET_TEST_NUM 38
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc37_w, 32
    SET_TEST_NUM 39
    CHECK_CSRS_UNCHANGED

    /* Test 40-42: vmadd.vv SEW=64 max */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc40_vd
    vle64.v v8, (t1)
    la t1, tc40_s1
    vle64.v v20, (t1)
    la t1, tc40_s2
    vle64.v v16, (t1)
    la t1, tc40_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vmadd.vv v8, v20, v16
    SET_TEST_NUM 40
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc40_exp, 32
    SET_TEST_NUM 41
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc40_w, 32
    SET_TEST_NUM 42
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_vd:
    .byte 0x00, 0x00, 0x00, 0x00
tc1_s1:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_s2:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc4_vd:
    .byte 0x64, 0xc8, 0x2c, 0x90
tc4_s1:
    .byte 0x01, 0x01, 0x01, 0x01
tc4_s2:
    .byte 0x01, 0x02, 0x03, 0x04
tc4_exp:
    .byte 0x65, 0xca, 0x2f, 0x94
tc4_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc7_vd:
    .byte 0x00, 0x00, 0x00, 0x00
tc7_s1:
    .byte 0xff, 0xff, 0x01, 0x00
tc7_s2:
    .byte 0x01, 0x01, 0xff, 0xff
tc7_exp:
    .byte 0x01, 0x01, 0xff, 0xff
tc7_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc10_vd:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc10_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc10_s2:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc10_exp:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc10_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc13_vd:
    .half 0x0064, 0x00c8, 0x012c, 0x0190
tc13_s1:
    .half 0x0001, 0x0001, 0x0001, 0x0001
tc13_s2:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc13_exp:
    .half 0x0065, 0x00ca, 0x012f, 0x0194
tc13_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc16_vd:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc16_s1:
    .half 0xffff, 0xffff, 0x0001, 0x0000
tc16_s2:
    .half 0x0001, 0x0001, 0xffff, 0xffff
tc16_exp:
    .half 0x0001, 0x0001, 0xffff, 0xffff
tc16_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 2
tc19_vd:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc19_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc19_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc19_exp:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc19_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc22_vd:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
tc22_s1:
    .word 0x00000001, 0x00000001, 0x00000001, 0x00000001
tc22_s2:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc22_exp:
    .word 0x00000065, 0x000000ca, 0x0000012f, 0x00000194
tc22_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc25_vd:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc25_s1:
    .word 0xffffffff, 0xffffffff, 0x00000001, 0x00000000
tc25_s2:
    .word 0x00000001, 0x00000001, 0xffffffff, 0xffffffff
tc25_exp:
    .word 0x00000001, 0x00000001, 0xffffffff, 0xffffffff
tc25_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc28_vd:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc28_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc28_exp:
    .word 0x00000014, 0x0000003c, 0x00000078, 0x000000c8
tc28_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 1
tc31_mask:
    .byte 10
.align 2
tc31_vd:
    .word 0x0000dead, 0x0000dead, 0x0000dead, 0x0000dead
tc31_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc31_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc31_exp:
    .word 0x0000dead, 0x0001bd6e, 0x0000dead, 0x00037adc
tc31_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc34_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc34_s1:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc34_s2:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc34_exp:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc34_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc37_vd:
    .dword 0x0000000000000064, 0x00000000000000c8, 0x000000000000012c, 0x0000000000000190
tc37_s1:
    .dword 0x0000000000000001, 0x0000000000000001, 0x0000000000000001, 0x0000000000000001
tc37_s2:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc37_exp:
    .dword 0x0000000000000065, 0x00000000000000ca, 0x000000000000012f, 0x0000000000000194
tc37_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc40_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc40_s1:
    .dword 0xffffffffffffffff, 0xffffffffffffffff, 0x0000000000000001, 0x0000000000000000
tc40_s2:
    .dword 0x0000000000000001, 0x0000000000000001, 0xffffffffffffffff, 0xffffffffffffffff
tc40_exp:
    .dword 0x0000000000000001, 0x0000000000000001, 0xffffffffffffffff, 0xffffffffffffffff
tc40_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

