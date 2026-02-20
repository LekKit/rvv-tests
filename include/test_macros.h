/* test_macros.h - Convenience macros for RVV tests
 *
 * These macros reduce boilerplate in test .S files.
 * Used by both hand-written and generated tests.
 */

#ifndef TEST_MACROS_H
#define TEST_MACROS_H

/* ============================================================
 * Test case management
 * ============================================================ */

/* Set current test number (used by fail handler) */
.macro SET_TEST_NUM num
    li s11, \num
.endm

/* Exit with failure: uses s11 as exit code */
.macro FAIL_TEST
    mv a0, s11
    li a7, 93
    ecall
.endm

/* Exit with explicit failure code */
.macro FAIL_CODE code
    li a0, \code
    li a7, 93
    ecall
.endm

/* Exit success */
.macro PASS_TEST
    li a0, 0
    li a7, 93
    ecall
.endm

/* ============================================================
 * Conditional failure macros
 * ============================================================ */

/* Fail if reg1 != reg2 */
.macro FAIL_IF_NE reg1, reg2
    beq \reg1, \reg2, 90192f
    FAIL_TEST
90192:
.endm

/* Fail if reg != 0 */
.macro FAIL_IF_NZ reg
    beqz \reg, 90193f
    FAIL_TEST
90193:
.endm

/* Fail if reg == 0 */
.macro FAIL_IF_Z reg
    bnez \reg, 90194f
    FAIL_TEST
90194:
.endm

/* ============================================================
 * Memory comparison macro
 * Compare nbytes at actual_label vs expected_label.
 * Fails current test (s11) on mismatch.
 * Clobbers: a0-a3, t0-t1, ra
 * ============================================================ */
.macro CHECK_MEM actual_label, expected_label, nbytes
    la   a1, \actual_label
    la   a2, \expected_label
    li   a3, \nbytes
    jal  ra, _mem_compare
    FAIL_IF_NZ a0
.endm

/* ============================================================
 * CSR check macros
 * ============================================================ */

/* Save V-extension CSRs into s-registers for later checking.
 * s0 = vstart, s1 = vxsat, s2 = vxrm, s3 = vl, s4 = vtype
 * Also s5 = fcsr (FP flags, since FP vector ops can modify fflags)
 */
.macro SAVE_CSRS
    csrr s0, vstart
    csrr s1, vxsat
    csrr s2, vxrm
    csrr s3, vl
    csrr s4, vtype
    csrr s5, fcsr
.endm

/* Check that vstart is 0 (all vector instructions reset vstart) */
.macro CHECK_VSTART_ZERO
    csrr t0, vstart
    FAIL_IF_NZ t0
.endm

/* Check vxsat unchanged */
.macro CHECK_VXSAT_UNCHANGED
    csrr t0, vxsat
    FAIL_IF_NE t0, s1
.endm

/* Check vxrm unchanged */
.macro CHECK_VXRM_UNCHANGED
    csrr t0, vxrm
    FAIL_IF_NE t0, s2
.endm

/* Check vl unchanged */
.macro CHECK_VL_UNCHANGED
    csrr t0, vl
    FAIL_IF_NE t0, s3
.endm

/* Check vtype unchanged */
.macro CHECK_VTYPE_UNCHANGED
    csrr t0, vtype
    FAIL_IF_NE t0, s4
.endm

/* Check fcsr unchanged (for non-FP operations) */
.macro CHECK_FCSR_UNCHANGED
    csrr t0, fcsr
    FAIL_IF_NE t0, s5
.endm

/* Full CSR side-effect check for non-FP, non-fixedpoint instructions */
.macro CHECK_CSRS_UNCHANGED
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_VL_UNCHANGED
    CHECK_VTYPE_UNCHANGED
    CHECK_FCSR_UNCHANGED
.endm

/* CSR check for FP instructions (fflags may change, that's expected) */
.macro CHECK_CSRS_UNCHANGED_FP
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_VL_UNCHANGED
    CHECK_VTYPE_UNCHANGED
    /* fcsr/fflags NOT checked - FP ops may legitimately modify it */
.endm

/* CSR check for fixed-point instructions (vxsat may change) */
.macro CHECK_CSRS_UNCHANGED_FIXEDPOINT
    CHECK_VSTART_ZERO
    /* vxsat NOT checked - fixed-point saturating ops may set it */
    CHECK_VXRM_UNCHANGED
    CHECK_VL_UNCHANGED
    CHECK_VTYPE_UNCHANGED
    CHECK_FCSR_UNCHANGED
.endm

/* ============================================================
 * Vector register witness macros
 * Initialize and verify "witness" vector registers that should
 * NOT be modified by the instruction under test.
 * ============================================================ */

/* Initialize a vector register with a known byte pattern.
 * The pattern is: every element = (pattern_byte repeated).
 * Requires vtype already set.
 * Clobbers: t0
 */
.macro INIT_WITNESS vreg, pattern_byte
    li t0, \pattern_byte
    vmv.v.x \vreg, t0
.endm

/* Check a witness register by storing to buf and comparing.
 * store_insn = the store instruction to use (e.g., vse32.v)
 * buf_label  = label of scratch buffer
 * exp_label  = label of expected data
 * nbytes     = number of bytes to compare
 * Clobbers: a0-a3, t0-t1, t4, ra
 */
.macro CHECK_WITNESS vreg, store_insn, buf_label, exp_label, nbytes
    la t4, \buf_label
    \store_insn \vreg, (t4)
    CHECK_MEM \buf_label, \exp_label, \nbytes
.endm

#endif /* TEST_MACROS_H */
