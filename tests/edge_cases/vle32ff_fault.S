/* Auto-generated test for vle32ff.v
 * Fault-only-first load: happy path + mmap boundary fault
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = happy ff: wrong data loaded
 *     2 = happy ff: vl changed (should stay 4)
 *     3 = mmap failed
 *     4 = munmap second page failed
 *     5 = ff vl == 0 (spec requires >= 1)
 *     6 = ff vl > 3 (loaded from unmapped page?)
 *     7 = ff vstart != 0 after vle32ff
 *     8 = ff element 0 wrong
 *     9 = ff element 1 wrong
 *    10 = ff element 2 wrong
 *    11 = munmap cleanup failed
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* ============================================================ */
    /* Sub-test A: Happy path -- vle32ff.v from fully mapped memory */
    /* ============================================================ */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_src
    vle32ff.v v8, (t1)
    SET_TEST_NUM 2
    csrr t2, vl
    li t3, 4
    FAIL_IF_NE t2, t3
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_src, 16

    /* ============================================================ */
    /* Sub-test B: Boundary fault -- mmap + munmap trick */
    /* ============================================================ */
    SET_TEST_NUM 3
    li a0, 0
    li a1, 8192
    li a2, 3
    li a3, 34
    li a4, -1
    li a5, 0
    li a7, 222
    ecall
    bltz a0, 90100f
    j 90101f
90100:
    FAIL_TEST
90101:
    mv s6, a0
    SET_TEST_NUM 4
    li t0, 4096
    add a0, s6, t0
    li a1, 4096
    li a7, 215
    ecall
    FAIL_IF_NZ a0
    li t0, 4084
    add s7, s6, t0
    li t1, 0xAAAA0001
    sw t1, 0(s7)
    li t1, 0xBBBB0002
    sw t1, 4(s7)
    li t1, 0xCCCC0003
    sw t1, 8(s7)
    vsetvli t0, x0, e32, m1, tu, mu
    li t1, 0xdeadbeef
    vmv.v.x v8, t1
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    vle32ff.v v8, (s7)
    csrr s8, vl
    SET_TEST_NUM 5
    beqz s8, 90102f
    j 90103f
90102:
    FAIL_TEST
90103:
    SET_TEST_NUM 6
    li t0, 3
    bgt s8, t0, 90104f
    j 90105f
90104:
    FAIL_TEST
90105:
    SET_TEST_NUM 7
    csrr t2, vstart
    FAIL_IF_NZ t2
    la t1, ff_resbuf
    vs1r.v v8, (t1)
    SET_TEST_NUM 8
    la a1, ff_resbuf
    la a2, ff_exp
    li a3, 4
    jal ra, _mem_compare
    FAIL_IF_NZ a0
    SET_TEST_NUM 9
    li t0, 2
    blt s8, t0, 90106f
    la a1, ff_resbuf
    addi a1, a1, 4
    la a2, ff_exp
    addi a2, a2, 4
    li a3, 4
    jal ra, _mem_compare
    FAIL_IF_NZ a0
90106:
    SET_TEST_NUM 10
    li t0, 3
    blt s8, t0, 90107f
    la a1, ff_resbuf
    addi a1, a1, 8
    la a2, ff_exp
    addi a2, a2, 8
    li a3, 4
    jal ra, _mem_compare
    FAIL_IF_NZ a0
90107:
    SET_TEST_NUM 11
    mv a0, s6
    li a1, 4096
    li a7, 215
    ecall
    FAIL_IF_NZ a0

    PASS_TEST

.data
.align 2
tc1_src:
    .word 0xaaaa0001, 0xbbbb0002, 0xcccc0003, 0xdddd0004
.align 2
ff_exp:
    .word 0xaaaa0001, 0xbbbb0002, 0xcccc0003
.align 4
ff_resbuf:
    .space 512

.align 4
result_buf:  .space 256
witness_buf: .space 256

