/* Auto-generated test for vfredosum.vs
 * FP reduction vfredosum.vs
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfredosum.vs e32 basic: wrong result
 *     2 = vfredosum.vs e32 basic: witness changed
 *     3 = vfredosum.vs e32 basic: CSR side-effect
 *     4 = vfredosum.vs e64 basic: wrong result
 *     5 = vfredosum.vs e64 basic: witness changed
 *     6 = vfredosum.vs e64 basic: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vfredosum.vs SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_vec
    vle32.v v16, (t1)
    la t1, tc1_init
    vle32.v v20, (t1)
    la t1, tc1_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfredosum.vs v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    SET_TEST_NUM 2
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 16
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED_FP

    /* Test 4-6: vfredosum.vs SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc4_vec
    vle64.v v16, (t1)
    la t1, tc4_init
    vle64.v v20, (t1)
    la t1, tc4_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vfredosum.vs v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 8
    SET_TEST_NUM 5
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 32
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_vec:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tc1_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc1_exp:
    .word 0x41200000, 0x00000000, 0x00000000, 0x00000000
tc1_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc4_vec:
    .dword 0x3ff0000000000000, 0x4000000000000000, 0x4008000000000000, 0x4010000000000000
tc4_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc4_exp:
    .dword 0x4024000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc4_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

