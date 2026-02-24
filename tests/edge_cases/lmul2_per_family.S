/* Auto-generated test for LMUL=2 per family
 * Tests register group boundary crossing (e32/m2) for vadd, vmul, vand, vsll, vmin across 16 elements (v8+v9)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vadd.vv e32/m2: boundary crossing
 *     2 = vmul.vv e32/m2: boundary crossing
 *     3 = vand.vv e32/m2: boundary crossing
 *     4 = vsll.vv e32/m2: boundary crossing
 *     5 = vminu.vv e32/m2: boundary crossing
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* LMUL=2 test: vadd.vv (e32/m2, vl=VLMAX) */
    vsetvli t0, x0, e32, m2, ta, ma
    mv s6, t0
    la t1, lm21_s2
    vle32.v v16, (t1)
    la t1, lm21_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vadd.vv v8, v16, v20
    la t1, lm21_res
    vse32.v v8, (t1)
    SET_TEST_NUM 1
    slli a3, s6, 2
    la a1, lm21_res
    la a2, lm21_exp
    call _mem_compare
    FAIL_IF_NZ a0

    /* LMUL=2 test: vmul.vv (e32/m2, vl=VLMAX) */
    vsetvli t0, x0, e32, m2, ta, ma
    mv s6, t0
    la t1, lm22_s2
    vle32.v v16, (t1)
    la t1, lm22_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmul.vv v8, v16, v20
    la t1, lm22_res
    vse32.v v8, (t1)
    SET_TEST_NUM 2
    slli a3, s6, 2
    la a1, lm22_res
    la a2, lm22_exp
    call _mem_compare
    FAIL_IF_NZ a0

    /* LMUL=2 test: vand.vv (e32/m2, vl=VLMAX) */
    vsetvli t0, x0, e32, m2, ta, ma
    mv s6, t0
    la t1, lm23_s2
    vle32.v v16, (t1)
    la t1, lm23_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vand.vv v8, v16, v20
    la t1, lm23_res
    vse32.v v8, (t1)
    SET_TEST_NUM 3
    slli a3, s6, 2
    la a1, lm23_res
    la a2, lm23_exp
    call _mem_compare
    FAIL_IF_NZ a0

    /* LMUL=2 test: vsll.vv (e32/m2, vl=VLMAX) */
    vsetvli t0, x0, e32, m2, ta, ma
    mv s6, t0
    la t1, lm24_s2
    vle32.v v16, (t1)
    la t1, lm24_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vsll.vv v8, v16, v20
    la t1, lm24_res
    vse32.v v8, (t1)
    SET_TEST_NUM 4
    slli a3, s6, 2
    la a1, lm24_res
    la a2, lm24_exp
    call _mem_compare
    FAIL_IF_NZ a0

    /* LMUL=2 test: vminu.vv (e32/m2, vl=VLMAX) */
    vsetvli t0, x0, e32, m2, ta, ma
    mv s6, t0
    la t1, lm25_s2
    vle32.v v16, (t1)
    la t1, lm25_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vminu.vv v8, v16, v20
    la t1, lm25_res
    vse32.v v8, (t1)
    SET_TEST_NUM 5
    slli a3, s6, 2
    la a1, lm25_res
    la a2, lm25_exp
    call _mem_compare
    FAIL_IF_NZ a0

    PASS_TEST

.data
.align 4
lm21_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190, 0x000001f4, 0x00000258, 0x000002bc, 0x00000320, 0x00000384, 0x000003e8, 0x0000044c, 0x000004b0, 0x00000514, 0x00000578, 0x000005dc, 0x00000640
lm21_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004, 0x00000005, 0x00000006, 0x00000007, 0x00000008, 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010
lm21_exp:
    .word 0x00000065, 0x000000ca, 0x0000012f, 0x00000194, 0x000001f9, 0x0000025e, 0x000002c3, 0x00000328, 0x0000038d, 0x000003f2, 0x00000457, 0x000004bc, 0x00000521, 0x00000586, 0x000005eb, 0x00000650
lm21_res:
    .space 64
.align 4
lm22_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190, 0x000001f4, 0x00000258, 0x000002bc, 0x00000320, 0x00000384, 0x000003e8, 0x0000044c, 0x000004b0, 0x00000514, 0x00000578, 0x000005dc, 0x00000640
lm22_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004, 0x00000005, 0x00000006, 0x00000007, 0x00000008, 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010
lm22_exp:
    .word 0x00000064, 0x00000190, 0x00000384, 0x00000640, 0x000009c4, 0x00000e10, 0x00001324, 0x00001900, 0x00001fa4, 0x00002710, 0x00002f44, 0x00003840, 0x00004204, 0x00004c90, 0x000057e4, 0x00006400
lm22_res:
    .space 64
.align 4
lm23_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190, 0x000001f4, 0x00000258, 0x000002bc, 0x00000320, 0x00000384, 0x000003e8, 0x0000044c, 0x000004b0, 0x00000514, 0x00000578, 0x000005dc, 0x00000640
lm23_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004, 0x00000005, 0x00000006, 0x00000007, 0x00000008, 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010
lm23_exp:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000004, 0x00000000, 0x00000004, 0x00000000, 0x00000000, 0x00000008, 0x00000008, 0x00000000, 0x00000004, 0x00000008, 0x0000000c, 0x00000000
lm23_res:
    .space 64
.align 4
lm24_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190, 0x000001f4, 0x00000258, 0x000002bc, 0x00000320, 0x00000384, 0x000003e8, 0x0000044c, 0x000004b0, 0x00000514, 0x00000578, 0x000005dc, 0x00000640
lm24_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004, 0x00000005, 0x00000006, 0x00000007, 0x00000008, 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010
lm24_exp:
    .word 0x000000c8, 0x00000320, 0x00000960, 0x00001900, 0x00003e80, 0x00009600, 0x00015e00, 0x00032000, 0x00070800, 0x000fa000, 0x00226000, 0x004b0000, 0x00a28000, 0x015e0000, 0x02ee0000, 0x06400000
lm24_res:
    .space 64
.align 4
lm25_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190, 0x000001f4, 0x00000258, 0x000002bc, 0x00000320, 0x00000384, 0x000003e8, 0x0000044c, 0x000004b0, 0x00000514, 0x00000578, 0x000005dc, 0x00000640
lm25_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004, 0x00000005, 0x00000006, 0x00000007, 0x00000008, 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010
lm25_exp:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004, 0x00000005, 0x00000006, 0x00000007, 0x00000008, 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010
lm25_res:
    .space 64

.align 4
result_buf:  .space 256
witness_buf: .space 256

