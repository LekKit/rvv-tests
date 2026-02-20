/* Auto-generated test for vwmaccsu.vv
 * Widening macc vwmaccsu.vv
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vwmaccsu.vv e8 basic: result
 *     2 = vwmaccsu.vv e16 basic: result
 *     3 = vwmaccsu.vv e32 basic: result
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1: vwmaccsu.vv SEW=8 */
    li t0, 4
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, tc1_vd
    vle16.v v8, (t1)
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_s1
    vle8.v v20, (t1)
    la t1, tc1_s2
    vle8.v v16, (t1)
    SAVE_CSRS
    vwmaccsu.vv v8, v20, v16
    SET_TEST_NUM 1
    vsetvli t0, t0, e16, m2, tu, mu
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 8
    CHECK_VSTART_ZERO

    /* Test 2: vwmaccsu.vv SEW=16 */
    li t0, 4
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc2_vd
    vle32.v v8, (t1)
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc2_s1
    vle16.v v20, (t1)
    la t1, tc2_s2
    vle16.v v16, (t1)
    SAVE_CSRS
    vwmaccsu.vv v8, v20, v16
    SET_TEST_NUM 2
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 16
    CHECK_VSTART_ZERO

    /* Test 3: vwmaccsu.vv SEW=32 */
    li t0, 4
    vsetvli t0, t0, e64, m2, tu, mu
    la t1, tc3_vd
    vle64.v v8, (t1)
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc3_s1
    vle32.v v20, (t1)
    la t1, tc3_s2
    vle32.v v16, (t1)
    SAVE_CSRS
    vwmaccsu.vv v8, v20, v16
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
tc1_s1:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_s2:
    .byte 0x0a, 0x14, 0x1e, 0x28
tc1_exp:
    .half 0x000a, 0x0028, 0x005a, 0x00a0
.align 2
tc2_vd:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc2_s1:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc2_s2:
    .half 0x000a, 0x0014, 0x001e, 0x0028
tc2_exp:
    .word 0x0000000a, 0x00000028, 0x0000005a, 0x000000a0
.align 3
tc3_vd:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc3_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc3_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tc3_exp:
    .dword 0x000000000000000a, 0x0000000000000028, 0x000000000000005a, 0x00000000000000a0

.align 4
result_buf:  .space 256
witness_buf: .space 256

