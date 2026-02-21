/* Auto-generated test for vadd.vv
 * Mask agnostic (ma): tu/ma mask=0b0101, verify only unmasked elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = mask agnostic: element 0 wrong
 *     2 = mask agnostic: element 2 wrong
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Mask agnostic: vadd.vv e32 tu/ma mask=0b0101 */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, ma
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    la t1, tc1_mask
    vlm.v v0, (t1)
    vadd.vv v8, v16, v20, v0.t
    la t1, result_buf
    vse32.v v8, (t1)
    SET_TEST_NUM 1
    la a1, result_buf
    la a2, tc1_exp0
    li a3, 4
    jal ra, _mem_compare
    FAIL_IF_NZ a0
    SET_TEST_NUM 2
    la a1, result_buf
    addi a1, a1, 8
    la a2, tc1_exp2
    li a3, 4
    jal ra, _mem_compare
    FAIL_IF_NZ a0

    PASS_TEST

.data
.align 1
tc1_mask:
    .byte 5
.align 2
tc1_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
tc1_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc1_exp0:
    .word 0x0000006e
tc1_exp2:
    .word 0x0000014a

.align 4
result_buf:  .space 256
witness_buf: .space 256

