/* Auto-generated test for GhostWrite CVE-2024-44067 check
 * Verifies mew=1 (reserved) vector store/load instructions trap as illegal instruction, preventing the GhostWrite physical-memory write vulnerability (CVE-2024-44067, T-Head C910/C920)
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = clone failed for vse128.v check
 *     2 = vse128.v (GhostWrite, 0x10028027) did NOT trap — VULNERABLE to CVE-2024-44067!
 *     3 = clone failed for vle128.v check
 *     4 = vle128.v (0x12028407, mew=1 load) did NOT trap — reserved encoding accepted
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* GhostWrite vulnerability check (CVE-2024-44067) */
    /* T-Head XuanTie C910/C920: vse128.v writes physical memory */
    /* Standard RVV 1.0: mew=1 (bit 28) is reserved → must SIGILL */

    /* SIGILL check: vse128.v v0, 0(t0) — GhostWrite store (CVE-2024-44067) */
    /* Encoding 0x10028027 — must trap as illegal instruction */
    SET_TEST_NUM 1
    SYS_CLONE
    bltz a0, gw_store_clone_fail
    beqz a0, gw_store_child
    j gw_store_parent
gw_store_clone_fail:
    FAIL_TEST

gw_store_child:
    li t0, 4
    vsetvli t0, t0, e8, m1, ta, ma
    la t0, gw_scratch
    /* Execute 0x10028027 — should SIGILL */
    .4byte 0x10028027
    /* If we reach here, instruction did NOT trap — exit child 0 */
    SYS_EXIT 0

gw_store_parent:
    mv s6, a0
    la s7, gw_store_wstatus
    SYS_WAIT4 s6, s7

    SET_TEST_NUM 2
    la t0, gw_store_wstatus
    lw t0, 0(t0)
    andi t0, t0, 0x7f
    li t1, 4
    FAIL_IF_NE t0, t1

    /* SIGILL check: vle128.v v8, 0(t0) — mew=1 load */
    /* Encoding 0x12028407 — must trap as illegal instruction */
    SET_TEST_NUM 3
    SYS_CLONE
    bltz a0, gw_load_clone_fail
    beqz a0, gw_load_child
    j gw_load_parent
gw_load_clone_fail:
    FAIL_TEST

gw_load_child:
    li t0, 4
    vsetvli t0, t0, e8, m1, ta, ma
    la t0, gw_scratch
    /* Execute 0x12028407 — should SIGILL */
    .4byte 0x12028407
    /* If we reach here, instruction did NOT trap — exit child 0 */
    SYS_EXIT 0

gw_load_parent:
    mv s6, a0
    la s7, gw_load_wstatus
    SYS_WAIT4 s6, s7

    SET_TEST_NUM 4
    la t0, gw_load_wstatus
    lw t0, 0(t0)
    andi t0, t0, 0x7f
    li t1, 4
    FAIL_IF_NE t0, t1

    PASS_TEST

.data
.align 4
gw_scratch:
    .space 64
.align 3
gw_store_wstatus:
    .word 0
gw_load_wstatus:
    .word 0

.align 4
result_buf:  .space 256
witness_buf: .space 256

