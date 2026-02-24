/* Auto-generated test for LMUL=4 and fractional LMUL
 * Tests vadd with e32/m4 (32 elements, 4 register groups) and e32/mf2 (fractional, reduced VLMAX)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vadd.vv e32/m4: boundary crossing
 *     2 = vadd.vv e32/mf2: fractional LMUL
 *     3 = vmul.vv e16/mf4: fractional LMUL
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* LMUL=4 vadd.vv (e32/m4, vl=VLMAX) */
    vsetvli t0, x0, e32, m4, ta, ma
    mv s6, t0
    la t1, lm4f1_s2
    vle32.v v16, (t1)
    la t1, lm4f1_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vadd.vv v8, v16, v20
    la t1, lm4f1_res
    vse32.v v8, (t1)
    SET_TEST_NUM 1
    slli a3, s6, 2
    la a1, lm4f1_res
    la a2, lm4f1_exp
    call _mem_compare
    FAIL_IF_NZ a0

    /* Fractional LMUL vadd.vv (e32/mf2) */
    vsetvli t0, x0, e32, mf2, ta, ma
    mv s6, t0
    la t1, lm4f2_s2
    vle32.v v16, (t1)
    la t1, lm4f2_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vadd.vv v8, v16, v20
    la t1, lm4f2_res
    vse32.v v8, (t1)
    SET_TEST_NUM 2
    slli a3, s6, 2
    la a1, lm4f2_res
    la a2, lm4f2_exp
    call _mem_compare
    FAIL_IF_NZ a0

    /* Fractional LMUL vmul.vv (e16/mf4) */
    vsetvli t0, x0, e16, mf4, ta, ma
    mv s6, t0
    la t1, lm4f3_s2
    vle16.v v16, (t1)
    la t1, lm4f3_s1
    vle16.v v20, (t1)
    SAVE_CSRS
    vmul.vv v8, v16, v20
    la t1, lm4f3_res
    vse16.v v8, (t1)
    SET_TEST_NUM 3
    slli a3, s6, 1
    la a1, lm4f3_res
    la a2, lm4f3_exp
    call _mem_compare
    FAIL_IF_NZ a0

    PASS_TEST

.data
.align 4
lm4f1_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028, 0x00000032, 0x0000003c, 0x00000046, 0x00000050, 0x0000005a, 0x00000064, 0x0000006e, 0x00000078, 0x00000082, 0x0000008c, 0x00000096, 0x000000a0, 0x000000aa, 0x000000b4, 0x000000be, 0x000000c8, 0x000000d2, 0x000000dc, 0x000000e6, 0x000000f0, 0x000000fa, 0x00000104, 0x0000010e, 0x00000118, 0x00000122, 0x0000012c, 0x00000136, 0x00000140
lm4f1_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004, 0x00000005, 0x00000006, 0x00000007, 0x00000008, 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c, 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010, 0x00000011, 0x00000012, 0x00000013, 0x00000014, 0x00000015, 0x00000016, 0x00000017, 0x00000018, 0x00000019, 0x0000001a, 0x0000001b, 0x0000001c, 0x0000001d, 0x0000001e, 0x0000001f, 0x00000020
lm4f1_exp:
    .word 0x0000000b, 0x00000016, 0x00000021, 0x0000002c, 0x00000037, 0x00000042, 0x0000004d, 0x00000058, 0x00000063, 0x0000006e, 0x00000079, 0x00000084, 0x0000008f, 0x0000009a, 0x000000a5, 0x000000b0, 0x000000bb, 0x000000c6, 0x000000d1, 0x000000dc, 0x000000e7, 0x000000f2, 0x000000fd, 0x00000108, 0x00000113, 0x0000011e, 0x00000129, 0x00000134, 0x0000013f, 0x0000014a, 0x00000155, 0x00000160
lm4f1_res:
    .space 128
.align 4
lm4f2_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
lm4f2_s1:
    .word 0x00000005, 0x0000000a, 0x0000000f, 0x00000014
lm4f2_exp:
    .word 0x00000069, 0x000000d2, 0x0000013b, 0x000001a4
lm4f2_res:
    .space 16
.align 2
lm4f3_s2:
    .half 0x0064, 0x00c8, 0x012c, 0x0190
lm4f3_s1:
    .half 0x0002, 0x0003, 0x0004, 0x0005
lm4f3_exp:
    .half 0x00c8, 0x0258, 0x04b0, 0x07d0
lm4f3_res:
    .space 8

.align 4
result_buf:  .space 256
witness_buf: .space 256

