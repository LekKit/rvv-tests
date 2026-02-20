/* Auto-generated test for vid.v
 * Vector element index
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vid.v e8: result
 *     2 = vid.v e16: result
 *     3 = vid.v e32: result
 *     4 = vid.v e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    SAVE_CSRS
    vid.v v8
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    SAVE_CSRS
    vid.v v8
    SET_TEST_NUM 2
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vid.v v8
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    SAVE_CSRS
    vid.v v8
    SET_TEST_NUM 4
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_exp:
    .byte 0x00, 0x01, 0x02, 0x03
.align 1
tc2_exp:
    .half 0x0000, 0x0001, 0x0002, 0x0003
.align 2
tc3_exp:
    .word 0x00000000, 0x00000001, 0x00000002, 0x00000003
.align 3
tc4_exp:
    .dword 0x0000000000000000, 0x0000000000000001, 0x0000000000000002, 0x0000000000000003

.align 4
result_buf:  .space 256
witness_buf: .space 256

