/* Auto-generated test for vstart CSR + vadd.vv
 * Verifies csrw/csrr vstart works and tests vadd.vv behavior with nonzero vstart (both SIGILL and skip-elements are valid per RVV 1.0 §3.7)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = csrw vstart, 2 caused trap (unexpected)
 *     2 = csrr vstart did not read back 2
 *     3 = clone() failed
 *     4 = vadd.vv with vstart=2: child exited 0 but result wrong (neither SIGILL nor correct skip-elements behavior)
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vstart nonzero test */
    /* RVV 1.0 §3.7: implementations may trap on nonzero vstart */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, vstart_vd
    vle32.v v8, (t1)
    la t1, vstart_s2
    vle32.v v16, (t1)
    la t1, vstart_s1
    vle32.v v20, (t1)
    SET_TEST_NUM 1
    li t0, 2
    csrw vstart, t0
    /* If we reach here, csrw vstart succeeded */
    SET_TEST_NUM 2
    csrr t1, vstart
    li t2, 2
    FAIL_IF_NE t1, t2
    csrw vstart, zero
    SET_TEST_NUM 3
    SYS_CLONE
    bltz a0, vstart_clone_fail
    beqz a0, vstart_child
    j vstart_parent
vstart_clone_fail:
    FAIL_TEST

vstart_child:
    /* Child: reload registers (not preserved across fork), */
    /* set vstart=2, execute vadd.vv */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, vstart_vd
    vle32.v v8, (t1)
    la t1, vstart_s2
    vle32.v v16, (t1)
    la t1, vstart_s1
    vle32.v v20, (t1)
    li t0, 2
    csrw vstart, t0
    /* This vadd.vv may SIGILL (legal) or skip elements 0-1 (also legal) */
    vadd.vv v8, v16, v20

    /* If we reach here, vadd.vv succeeded — store result and exit 0 */
    la t1, vstart_result
    vse32.v v8, (t1)
    SYS_EXIT 0

vstart_parent:
    mv s6, a0
    la s7, vstart_wstatus
    SYS_WAIT4 s6, s7

    /* Check child outcome */
    la t0, vstart_wstatus
    lw t0, 0(t0)

    /* WIFEXITED: bits [15:8] = exit status, bits [6:0] = 0 */
    /* WIFSIGNALED: bits [6:0] = signal number */
    andi t1, t0, 0x7f
    beqz t1, vstart_exited

    /* Child was signaled — check it was SIGILL (4) */
    li t2, 4
    beq t1, t2, vstart_sigill_ok
    /* Killed by unexpected signal — fail */
    SET_TEST_NUM 4
    FAIL_TEST

vstart_sigill_ok:
    /* SIGILL from vadd.vv with vstart=2 — legal per spec, pass */
    j vstart_done

vstart_exited:
    /* Child exited normally — vadd.vv ran, verify result */
    SET_TEST_NUM 4
    CHECK_MEM vstart_result, vstart_exp, 16

vstart_done:

    PASS_TEST

.data
.align 2
vstart_vd:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
vstart_s2:
    .word 0x00000064, 0x000000c8, 0x0000012c, 0x00000190
vstart_s1:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
vstart_exp:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0x0000014a, 0x000001b8
.align 3
vstart_wstatus:
    .word 0
vstart_result:
    .space 16

.align 4
result_buf:  .space 256
witness_buf: .space 256

