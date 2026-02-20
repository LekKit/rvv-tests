/* Auto-generated test for vfredmax.vs
 * FP reduction vfredmax.vs
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfredmax.vs e32 basic: wrong result
 *     2 = vfredmax.vs e32 basic: CSR side-effect
 *     3 = vfredmax.vs e64 basic: wrong result
 *     4 = vfredmax.vs e64 basic: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vfredmax.vs SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_vec
    vle32.v v16, (t1)
    la t1, tc1_init
    vle32.v v20, (t1)
    SAVE_CSRS
    vfredmax.vs v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED_FP

    /* Test 3-4: vfredmax.vs SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc3_vec
    vle64.v v16, (t1)
    la t1, tc3_init
    vle64.v v20, (t1)
    SAVE_CSRS
    vfredmax.vs v8, v16, v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 8
    SET_TEST_NUM 4
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_vec:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc1_exp:
    .word 0x40800000, 0x00000000, 0x00000000, 0x00000000
.align 3
tc3_vec:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc3_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc3_exp:
    .dword 0x4010000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

