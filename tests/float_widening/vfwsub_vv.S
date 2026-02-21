/* Auto-generated test for vfwsub.vv
 * Widening FP vfwsub.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfwsub.vv basic: result
 *     2 = vfwsub.vv basic: CSR side-effect
 *     3 = vfwsub.vv negpos: result
 *     4 = vfwsub.vv negpos: CSR side-effect
 *     5 = vfwsub.vv e32→e64 basic masked: wrong result
 *     6 = vfwsub.vv e32→e64 basic masked: witness changed
 *     7 = vfwsub.vv e32→e64 basic masked: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vfwsub.vv v8, v16, v20
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED_FP
    SET_TEST_NUM 1
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_s2
    vle32.v v16, (t1)
    la t1, tc3_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vfwsub.vv v8, v16, v20
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED_FP
    SET_TEST_NUM 3
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 32

    /* Test 5-7: vfwsub.vv SEW=32 basic (masked) */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc5_vd
    vle64.v v8, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc5_s2
    vle32.v v16, (t1)
    la t1, tc5_s1
    vle32.v v20, (t1)
    la t1, tc5_mask
    vlm.v v0, (t1)
    la t1, tc5_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfwsub.vv v8, v16, v20, v0.t
    SET_TEST_NUM 6
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc5_w, 16
    SET_TEST_NUM 7
    CHECK_CSRS_UNCHANGED_FP
    SET_TEST_NUM 5
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 32

    PASS_TEST

.data
.align 3
tc1_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_s1:
    .word 0x3f000000, 0x3e800000, 0x3e000000, 0x3d800000
tc1_exp:
    .dword 0x3fe0000000000000, 0x3ffc000000000000, 0x4007000000000000, 0x400f800000000000
.align 3
tc3_s2:
    .word 0xbf800000, 0xc0000000, 0x40400000, 0x40800000
tc3_s1:
    .word 0x3f800000, 0x40000000, 0xc0400000, 0xc0800000
tc3_exp:
    .dword 0xc000000000000000, 0xc010000000000000, 0x4018000000000000, 0x4020000000000000
.align 1
tc5_mask:
    .byte 10
.align 3
tc5_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc5_s1:
    .word 0x3f000000, 0x3e800000, 0x3e000000, 0x3d800000
tc5_vd:
    .dword 0x408f380000000000, 0x408f380000000000, 0x408f380000000000, 0x408f380000000000
tc5_exp:
    .dword 0x408f380000000000, 0x3ffc000000000000, 0x408f380000000000, 0x400f800000000000
tc5_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0

.align 4
result_buf:  .space 256
witness_buf: .space 256

