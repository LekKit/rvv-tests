/* Auto-generated test for vmflt.vv
 * FP compare vmflt.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmflt.vv e32 basic: wrong mask result
 *     2 = vmflt.vv e32 basic: CSR side-effect
 *     3 = vmflt.vv e32 negpos: wrong mask result
 *     4 = vmflt.vv e32 negpos: CSR side-effect
 *     5 = vmflt.vv e32 zeros: wrong mask result
 *     6 = vmflt.vv e32 zeros: CSR side-effect
 *     7 = vmflt.vv e64 basic: wrong mask result
 *     8 = vmflt.vv e64 basic: CSR side-effect
 *     9 = vmflt.vv e64 negpos: wrong mask result
 *    10 = vmflt.vv e64 negpos: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vmflt.vv SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmflt.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    /* Test 3-4: vmflt.vv SEW=32 negpos */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_s2
    vle32.v v16, (t1)
    la t1, tc3_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmflt.vv v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED

    /* Test 5-6: vmflt.vv SEW=32 zeros */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc5_s2
    vle32.v v16, (t1)
    la t1, tc5_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmflt.vv v8, v16, v20
    SET_TEST_NUM 5
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-8: vmflt.vv SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc7_s2
    vle64.v v16, (t1)
    la t1, tc7_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vmflt.vv v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 0
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 8
    CHECK_CSRS_UNCHANGED

    /* Test 9-10: vmflt.vv SEW=64 negpos */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc9_s2
    vle64.v v16, (t1)
    la t1, tc9_s1
    vle64.v v20, (t1)
    SAVE_CSRS
    vmflt.vv v8, v16, v20
    SET_TEST_NUM 9
    la t1, result_buf
    vsm.v v8, (t1)
    lbu t2, 0(t1)
    andi t2, t2, 15
    li t3, 3
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 10
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_s1:
    .word 0x3f000000, 0x3e800000, 0x3e000000, 0x3d800000
.align 2
tc3_s2:
    .word 0xbf800000, 0xc0000000, 0x40400000, 0x40800000
tc3_s1:
    .word 0x3f800000, 0x40000000, 0xc0400000, 0xc0800000
.align 2
tc5_s2:
    .word 0x00000000, 0x80000000, 0x00000000, 0x80000000
tc5_s1:
    .word 0x00000000, 0x00000000, 0x80000000, 0x80000000
.align 3
tc7_s2:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc7_s1:
    .dword 0x3fe0000000000000, 0x3fd0000000000000, 0x3fc0000000000000, 0x3fb0000000000000
.align 3
tc9_s2:
    .dword 0xbff0000000000000, 0xc000000000000000, 0x4008000000000000, 0x4010000000000000
tc9_s1:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0xc008000000000000, 0xc010000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

