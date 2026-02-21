/* Auto-generated test for vsoxei32.v
 * Ordered indexed store: duplicate indices, later element overwrites earlier at same offset
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = scatter ordered: wrong memory content
 *     2 = scatter ordered: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Ordered indexed store with duplicate indices */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_data
    vle32.v v8, (t1)
    la t1, tc1_idx
    vle32.v v16, (t1)
    SAVE_CSRS
    la t1, tc1_target
    vsoxei32.v v8, (t1), v16
    SET_TEST_NUM 1
    CHECK_MEM tc1_target, tc1_exp, 16
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_data:
    .word 0x11111111, 0x22222222, 0x33333333, 0x44444444
tc1_idx:
    .word 0x00000008, 0x00000004, 0x00000004, 0x00000000
tc1_target:
    .word 0xaaaaaaaa, 0xaaaaaaaa, 0xaaaaaaaa, 0xaaaaaaaa
tc1_exp:
    .word 0x44444444, 0x33333333, 0x11111111, 0xaaaaaaaa

.align 4
result_buf:  .space 256
witness_buf: .space 256

