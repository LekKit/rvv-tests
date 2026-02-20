/* Auto-generated test for vwmacc.vx
 * Widening macc vwmacc.vx
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwmacc.vx e8 basic: result
 *     2 = vwmacc.vx e16 basic: result
 *     3 = vwmacc.vx e32 basic: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vwmacc.vx SEW=8 */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc1_vd
    vle16.v v8, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s2
    vle8.v v16, (t1)
    li a0, 0x02
    SAVE_CSRS
    vwmacc.vx v8, a0, v16
    SET_TEST_NUM 1
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8
    CHECK_VSTART_ZERO

    /* Test 2: vwmacc.vx SEW=16 */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc2_vd
    vle32.v v8, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_s2
    vle16.v v16, (t1)
    li a0, 0x0002
    SAVE_CSRS
    vwmacc.vx v8, a0, v16
    SET_TEST_NUM 2
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 16
    CHECK_VSTART_ZERO

    /* Test 3: vwmacc.vx SEW=32 */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc3_vd
    vle64.v v8, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_s2
    vle32.v v16, (t1)
    li a0, 0x00000002
    SAVE_CSRS
    vwmacc.vx v8, a0, v16
    SET_TEST_NUM 3
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 32
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 1
tc1_vd:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc1_s2:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .half 0x0014, 0x0028, 0x003c, 0x0050
.align 2
tc2_vd:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc2_s2:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc2_exp:
    .word 0x00000014, 0x00000028, 0x0000003c, 0x00000050
.align 3
tc3_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc3_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc3_exp:
    .dword 0x0000000000000014, 0x0000000000000028, 0x000000000000003c, 0x0000000000000050

.align 4
result_buf:  .space 256
witness_buf: .space 256

