/* Auto-generated test for vsetvli / vadd.vv / vmv1r.v
 * vill bit: illegal vtype causes SIGILL on vtype-dependent instructions; vmv1r.v works after restoring valid vtype
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vill bit not set after illegal vtype
 *     2 = vl not zero after illegal vtype
 *     3 = clone() syscall failed
 *     4 = vadd.vv did not trap with illegal vtype (child not killed by SIGILL)
 *     5 = vmv1r.v data wrong after restore
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vill trap test: illegal vtype → SIGILL on vadd.vv */
    /* Uses fork/wait to verify trap without userspace signal handler */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, vill_src
    vle32.v v16, (t1)
    vsetvli t0, x0, e64, mf2, ta, ma
    SET_TEST_NUM 1
    csrr t0, vtype
    bgez t0, 90200f
    j 90201f
90200:
    FAIL_TEST
90201:
    SET_TEST_NUM 2
    csrr t0, vl
    FAIL_IF_NZ t0
    SET_TEST_NUM 3

    /* Fork child process */
    SYS_CLONE
    bltz a0, 90202f
    beqz a0, 90203f
    j 90204f
90202:
    FAIL_TEST

90203:
    /* Child: set illegal vtype and execute vadd.vv */
    /* (must set vtype ourselves — kernel does not preserve */
    /* vtype/vl across fork on all platforms) */
    vsetvli t0, x0, e64, mf2, ta, ma

    /* This vadd.vv should SIGILL → kernel kills child */
    vadd.vv v8, v16, v20

    /* If we reach here, vadd.vv did not trap — exit child */
    SYS_EXIT 0

90204:
    /* Parent: wait for child via wait4() */
    mv s6, a0
    la s7, vill_wstatus
    SYS_WAIT4 s6, s7

    /* Check WTERMSIG(wstatus) == SIGILL (4) */
    SET_TEST_NUM 4
    la t0, vill_wstatus
    lw t0, 0(t0)
    andi t0, t0, 0x7f
    li t1, 4
    FAIL_IF_NE t0, t1

    /* Restore valid vtype and re-load v16 (vector regs may not */
    /* survive clone+wait due to lazy vector context management) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, vill_src
    vle32.v v16, (t1)
    vmv1r.v v8, v16
    SET_TEST_NUM 5
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, vill_src, 16

    PASS_TEST

.data
.align 2
vill_src:
    .word 0x11111111, 0x22222222, 0x33333333, 0x44444444
.align 3
vill_wstatus:
.word 0

.align 4
result_buf:  .space 256
witness_buf: .space 256

