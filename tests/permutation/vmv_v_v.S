/* Auto-generated test for vmv.v.v
 * Vector-vector move
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmv.v.v e8: result
 *     2 = vmv.v.v e16: result
 *     3 = vmv.v.v e32: result
 *     4 = vmv.v.v e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    vle8.v v16, (t1)
    SAVE_CSRS
    vmv.v.v v8, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_src, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_src
    vle16.v v16, (t1)
    SAVE_CSRS
    vmv.v.v v8, v16
    SET_TEST_NUM 2
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_src, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_src
    vle32.v v16, (t1)
    SAVE_CSRS
    vmv.v.v v8, v16
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc3_src, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc4_src
    vle64.v v16, (t1)
    SAVE_CSRS
    vmv.v.v v8, v16
    SET_TEST_NUM 4
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc4_src, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_src:
    .byte 0x0a, 0x14, 0x1e, 0x28
.align 1
tc2_src:
    .half 0x000a, 0x0014, 0x001e, 0x0028
.align 2
tc3_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
.align 3
tc4_src:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028

.align 4
result_buf:  .space 256
witness_buf: .space 256

