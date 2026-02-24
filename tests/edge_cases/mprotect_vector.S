/* Auto-generated test for mprotect + vector load/store
 * Tests that vector loads work on RO pages and vector stores to RO pages fault (SIGSEGV) after mprotect RW→RO.
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = mmap failed
 *     2 = mprotect failed
 *     3 = vector load from RO page: data mismatch
 *     4 = clone failed for RO store check
 *     5 = vector store to RO page did NOT fault
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* mprotect interaction with vector load/store */
    SET_TEST_NUM 1
    SYS_MMAP 0, 4096, 3, 34, -1, 0
    bltz a0, mprot_fail
    mv s6, a0
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    la t1, mprot_data
    vle32.v v16, (t1)
    vse32.v v16, (s6)
    SET_TEST_NUM 2
    SYS_MPROTECT s6, 4096, 1
    bltz a0, mprot_fail
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    vle32.v v8, (s6)
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, mprot_data, 16
    SET_TEST_NUM 4
    SYS_CLONE
    bltz a0, mprot_fail
    beqz a0, mprot_child
    j mprot_parent
mprot_child:
    /* Child: attempt vector store to RO page → should SIGSEGV */
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    vmv.v.i v16, 0
    vse32.v v16, (s6)
    SYS_EXIT 0

mprot_parent:
    mv s7, a0
    la s8, mprot_wstatus
    SYS_WAIT4 s7, s8
    SET_TEST_NUM 5
    la t0, mprot_wstatus
    lw t0, 0(t0)
    andi t0, t0, 0x7f
    li t1, 11
    FAIL_IF_NE t0, t1
    SYS_MPROTECT s6, 4096, 3
    SYS_MUNMAP s6, 4096
    j mprot_done
mprot_fail:
    FAIL_TEST
mprot_done:

    PASS_TEST

.data
.align 4
mprot_data:
    .word 0xaaaa1111, 0xbbbb2222, 0xcccc3333, 0xdddd4444
.align 3
mprot_wstatus:
    .word 0

.align 4
result_buf:  .space 256
witness_buf: .space 256

