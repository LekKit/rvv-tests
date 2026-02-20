/* Auto-generated test for vfcvt.f.x.v
 * FP conversion vfcvt.f.x.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vfcvt.f.x.v e32: result
 *     2 = vfcvt.f.x.v e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    vle32.v v16, (t1)
    SAVE_CSRS
    vfcvt.f.x.v v8, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 16
    CHECK_CSRS_UNCHANGED_FP

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc2_src
    vle64.v v16, (t1)
    SAVE_CSRS
    vfcvt.f.x.v v8, v16
    SET_TEST_NUM 2
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 32
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0x00000000, 0x00000001, 0x0000002a, 0x7fffffff
tc1_exp:
    .word 0x00000000, 0x3f800000, 0x42280000, 0x4f000000
.align 3
tc2_src:
    .dword 0x0000000000000000, 0x0000000000000001, 0x000000000000002a, 0x7fffffffffffffff
tc2_exp:
    .dword 0x0000000000000000, 0x3ff0000000000000, 0x4045000000000000, 0x43e0000000000000

.align 4
result_buf:  .space 256
witness_buf: .space 256

