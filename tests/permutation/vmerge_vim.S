/* Auto-generated test for vmerge.vim
 * Vector merge with immediate and mask
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmerge.vim e8: result
 *     2 = vmerge.vim e16: result
 *     3 = vmerge.vim e32: result
 *     4 = vmerge.vim e64: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    la t1, tc1_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmerge.vim v8, v16, 5, v0
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 4
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_s2
    vle16.v v16, (t1)
    la t1, tc2_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmerge.vim v8, v16, 5, v0
    SET_TEST_NUM 2
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 8
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_s2
    vle32.v v16, (t1)
    la t1, tc3_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmerge.vim v8, v16, 5, v0
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 16
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc4_s2
    vle64.v v16, (t1)
    la t1, tc4_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    vmerge.vim v8, v16, 5, v0
    SET_TEST_NUM 4
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_mask:
    .byte 5
.align 1
tc1_s2:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .byte 0x05, 0x14, 0x05, 0x28
.align 1
tc2_mask:
    .byte 5
.align 1
tc2_s2:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc2_exp:
    .half 0x0005, 0x0014, 0x0005, 0x0028
.align 1
tc3_mask:
    .byte 5
.align 2
tc3_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc3_exp:
    .word 0x00000005, 0x00000014, 0x00000005, 0x00000028
.align 1
tc4_mask:
    .byte 5
.align 3
tc4_s2:
    .dword 0x000000000000000a, 0x0000000000000014, 0x000000000000001e, 0x0000000000000028
tc4_exp:
    .dword 0x0000000000000005, 0x0000000000000014, 0x0000000000000005, 0x0000000000000028

.align 4
result_buf:  .space 256
witness_buf: .space 256

