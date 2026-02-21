/* riscv_test.h - Test framework for RVV 1.0 instruction tests
 *
 * Each .S test file includes this header. It provides:
 *   - _start entry point (sets up stack, jumps to _test_start)
 *   - _mem_compare utility function
 *
 * Convention:
 *   s11 = current test case number (set before each check)
 *   Exit code 0 = all tests passed
 *   Exit code N = test case N failed
 *
 * Register usage for test infrastructure:
 *   s0-s10  : saved CSR snapshots and temporaries across test cases
 *   s11     : current test number
 *   a0-a5   : arguments for utility functions
 *   t0-t6   : temporaries within test cases
 *   ra      : return address for utility calls
 *   sp      : stack pointer (set up at entry)
 *
 * Scalar registers (x1-x31) are assumed correct (GC extension).
 * Vector state is what we are testing.
 */

#ifndef RISCV_TEST_H
#define RISCV_TEST_H

/* ---- Entry point & utility code ---- */
.text
.globl _start
.align 2
_start:
    /* Set up a small stack (we need sp for function calls) */
    la sp, _stack_top
    j _test_start

/* ---- Memory compare utility ----
 * a1 = pointer to actual data
 * a2 = pointer to expected data
 * a3 = number of bytes to compare
 * Returns: a0 = 0 if equal, a0 = 1 if mismatch
 * Clobbers: t0, t1, a1, a2, a3
 */
.align 2
_mem_compare:
    beqz a3, .Lmem_cmp_eq
    lbu  t0, 0(a1)
    lbu  t1, 0(a2)
    bne  t0, t1, .Lmem_cmp_ne
    addi a1, a1, 1
    addi a2, a2, 1
    addi a3, a3, -1
    j    _mem_compare
.Lmem_cmp_eq:
    li a0, 0
    ret
.Lmem_cmp_ne:
    li a0, 1
    ret

.align 2
_test_start:

/* ---- Stack space (4096 bytes) ---- */
.pushsection .bss
.align 4
_stack_bottom:
    .space 4096
_stack_top:
.popsection

#endif /* RISCV_TEST_H */
