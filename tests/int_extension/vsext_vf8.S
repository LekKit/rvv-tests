/* Auto-generated test for vsext.vf8
 * Integer extension vsext.vf8
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsext.vf8 e8→e64 basic: result
 *     2 = vsext.vf8 e8→e64 basic: CSR side-effect
 *     3 = vsext.vf8 e8→e64 max: result
 *     4 = vsext.vf8 e8→e64 max: CSR side-effect
 *     5 = vsext.vf8 e8→e64 high: result
 *     6 = vsext.vf8 e8→e64 high: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-2: vsext.vf8 e8→e64 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e64, m8, tu, mu
    SAVE_CSRS
    vsext.vf8 v8, v16
    SET_TEST_NUM 1
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32
    SET_TEST_NUM 2
    CHECK_VSTART_ZERO

    /* Test 3-4: vsext.vf8 e8→e64 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc3_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e64, m8, tu, mu
    SAVE_CSRS
    vsext.vf8 v8, v16
    SET_TEST_NUM 3
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 32
    SET_TEST_NUM 4
    CHECK_VSTART_ZERO

    /* Test 5-6: vsext.vf8 e8→e64 high */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc5_src
    vle8.v v16, (t1)
    vsetvli t0, t0, e64, m8, tu, mu
    SAVE_CSRS
    vsext.vf8 v8, v16
    SET_TEST_NUM 5
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc5_exp, 32
    SET_TEST_NUM 6
    CHECK_VSTART_ZERO

    PASS_TEST

.data
.align 3
tc1_src:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_exp:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
.align 3
tc3_src:
    .byte 0xff, 0xfe, 0x7f, 0x00
tc3_exp:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0x000000000000007f, 0x0000000000000000
.align 3
tc5_src:
    .byte 0xff, 0x80, 0x81, 0x01
tc5_exp:
    .dword 0xffffffffffffffff, 0xffffffffffffff80, 0xffffffffffffff81, 0x0000000000000001

.align 4
result_buf:  .space 256
witness_buf: .space 256

