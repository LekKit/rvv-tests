/* Auto-generated test for vfsub.vv
 * FP vfsub.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfsub.vv e32 basic: wrong result
 *     2 = vfsub.vv e32 basic: witness changed
 *     3 = vfsub.vv e32 basic: CSR side-effect
 *     4 = vfsub.vv e32 negpos: wrong result
 *     5 = vfsub.vv e32 negpos: witness changed
 *     6 = vfsub.vv e32 negpos: CSR side-effect
 *     7 = vfsub.vv e32 zeros: wrong result
 *     8 = vfsub.vv e32 zeros: witness changed
 *     9 = vfsub.vv e32 zeros: CSR side-effect
 *    10 = vfsub.vv e32 basic overlap(vd=vs2): wrong result
 *    11 = vfsub.vv e32 basic overlap(vd=vs2): witness changed
 *    12 = vfsub.vv e32 basic overlap(vd=vs2): CSR side-effect
 *    13 = vfsub.vv e32 basic masked: wrong result
 *    14 = vfsub.vv e32 basic masked: witness changed
 *    15 = vfsub.vv e32 basic masked: CSR side-effect
 *    16 = vfsub.vv e64 basic: wrong result
 *    17 = vfsub.vv e64 basic: witness changed
 *    18 = vfsub.vv e64 basic: CSR side-effect
 *    19 = vfsub.vv e64 negpos: wrong result
 *    20 = vfsub.vv e64 negpos: witness changed
 *    21 = vfsub.vv e64 negpos: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vfsub.vv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    la t1, tc1_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfsub.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    SET_TEST_NUM 2
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 16
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED_FP

    /* Test 4-6: vfsub.vv SEW=32 negpos */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc4_s2
    vle32.v v16, (t1)
    la t1, tc4_s1
    vle32.v v20, (t1)
    la t1, tc4_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfsub.vv v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 16
    SET_TEST_NUM 5
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 16
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED_FP

    /* Test 7-9: vfsub.vv SEW=32 zeros */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc7_s2
    vle32.v v16, (t1)
    la t1, tc7_s1
    vle32.v v20, (t1)
    la t1, tc7_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfsub.vv v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 16
    SET_TEST_NUM 8
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 16
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED_FP

    /* Test 10-12: vfsub.vv SEW=32 basic (vd=vs2) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc10_s2
    vle32.v v8, (t1)
    la t1, tc10_s1
    vle32.v v20, (t1)
    la t1, tc10_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfsub.vv v8, v8, v20
    SET_TEST_NUM 10
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 16
    SET_TEST_NUM 11
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 16
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED_FP

    /* Test 13-15: vfsub.vv SEW=32 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc13_s2
    vle32.v v16, (t1)
    la t1, tc13_s1
    vle32.v v20, (t1)
    la t1, tc13_vd
    vle32.v v8, (t1)
    la t1, tc13_mask
    vlm.v v0, (t1)
    la t1, tc13_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfsub.vv v8, v16, v20, v0.t
    SET_TEST_NUM 13
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 16
    SET_TEST_NUM 14
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 16
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED_FP

    /* Test 16-18: vfsub.vv SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc16_s2
    vle64.v v16, (t1)
    la t1, tc16_s1
    vle64.v v20, (t1)
    la t1, tc16_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vfsub.vv v8, v16, v20
    SET_TEST_NUM 16
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 32
    SET_TEST_NUM 17
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc16_w, 32
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED_FP

    /* Test 19-21: vfsub.vv SEW=64 negpos */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc19_s2
    vle64.v v16, (t1)
    la t1, tc19_s1
    vle64.v v20, (t1)
    la t1, tc19_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vfsub.vv v8, v16, v20
    SET_TEST_NUM 19
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 32
    SET_TEST_NUM 20
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc19_w, 32
    SET_TEST_NUM 21
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_s1:
    .word 0x3f000000, 0x3e800000, 0x3e000000, 0x3d800000
tc1_exp:
    .word 0x3f000000, 0x3fe00000, 0x40380000, 0x407c0000
tc1_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc4_s2:
    .word 0xbf800000, 0xc0000000, 0x40400000, 0x40800000
tc4_s1:
    .word 0x3f800000, 0x40000000, 0xc0400000, 0xc0800000
tc4_exp:
    .word 0xc0000000, 0xc0800000, 0x40c00000, 0x41000000
tc4_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc7_s2:
    .word 0x00000000, 0x80000000, 0x00000000, 0x80000000
tc7_s1:
    .word 0x00000000, 0x00000000, 0x80000000, 0x80000000
tc7_exp:
    .word 0x00000000, 0x80000000, 0x00000000, 0x00000000
tc7_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc10_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc10_s1:
    .word 0x3f000000, 0x3e800000, 0x3e000000, 0x3d800000
tc10_exp:
    .word 0x3f000000, 0x3fe00000, 0x40380000, 0x407c0000
tc10_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 1
tc13_mask:
    .byte 10
.align 2
tc13_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc13_s1:
    .word 0x3f000000, 0x3e800000, 0x3e000000, 0x3d800000
tc13_vd:
    .word 0x42c80000, 0x42c80000, 0x42c80000, 0x42c80000
tc13_exp:
    .word 0x42c80000, 0x3fe00000, 0x42c80000, 0x407c0000
tc13_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc16_s2:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc16_s1:
    .dword 0x3fe0000000000000, 0x3fd0000000000000, 0x3fc0000000000000, 0x3fb0000000000000
tc16_exp:
    .dword 0x3fe0000000000000, 0x3ffc000000000000, 0x4007000000000000, 0x400f800000000000
tc16_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc19_s2:
    .dword 0xbff0000000000000, 0xc000000000000000, 0x4008000000000000, 0x4010000000000000
tc19_s1:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0xc008000000000000, 0xc010000000000000
tc19_exp:
    .dword 0xc000000000000000, 0xc010000000000000, 0x4018000000000000, 0x4020000000000000
tc19_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

