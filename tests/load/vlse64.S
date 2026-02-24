/* Auto-generated test for vlse64.v
 * Strided load of 64-bit elements
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vlse64.v stride=16: result
 *     2 = vlse64.v stride=16 masked: result
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_data
    li t2, 16
    SAVE_CSRS
    vlse64.v v8, (t1), t2
    SET_TEST_NUM 1
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32
    CHECK_CSRS_UNCHANGED

    /* Masked strided load: mask=0b0101 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc2_vd
    vle64.v v8, (t1)
    la t1, tc2_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    la t1, tc2_data
    li t2, 16
    vlse64.v v8, (t1), t2, v0.t
    SET_TEST_NUM 2
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_data:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004, 0x0000000000000005, 0x0000000000000006, 0x0000000000000007, 0x0000000000000008
tc1_exp:
    .dword 0x0000000000000001, 0x0000000000000003, 0x0000000000000005, 0x0000000000000007
.align 1
tc2_mask:
    .byte 5
.align 3
tc2_vd:
    .dword 0x00000000000000dd, 0x00000000000000dd, 0x00000000000000dd, 0x00000000000000dd
tc2_data:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004, 0x0000000000000005, 0x0000000000000006, 0x0000000000000007, 0x0000000000000008
tc2_exp:
    .dword 0x0000000000000001, 0x00000000000000dd, 0x0000000000000005, 0x00000000000000dd

.align 4
result_buf:  .space 256
witness_buf: .space 256

