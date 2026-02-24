/* Auto-generated test for Reserved RVV encoding checks
 * Verifies that structurally reserved RVV 1.0 load/store encodings trap as illegal instruction (SIGILL). Tests mew=1 (strided store and indexed load) and reserved lumop/sumop values.
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = clone failed for mew=1 strided store
 *     2 = mew=1 strided store (0x1A628027) did NOT trap
 *     3 = clone failed for reserved lumop
 *     4 = reserved lumop=00001 (0x02128407) did NOT trap
 *     5 = clone failed for mew=1 indexed load
 *     6 = mew=1 indexed load (0x17028407) did NOT trap
 *     7 = clone failed for reserved sumop
 *     8 = reserved sumop=00001 (0x02128427) did NOT trap
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Reserved RVV 1.0 load/store encoding checks */
    /* These test structurally reserved fields (mew, lumop, sumop), */
    /* NOT OP-V funct6 slots (which are claimed by newer extensions). */

    /* SIGILL check: mew=1 strided store */
    /* Encoding 0x1A628027 — must trap as illegal instruction */
    SET_TEST_NUM 1
    SYS_CLONE
    bltz a0, re_smew_clone_fail
    beqz a0, re_smew_child
    j re_smew_parent
re_smew_clone_fail:
    FAIL_TEST

re_smew_child:
    li t0, 4
    vsetvli t0, t0, e8, m1, ta, ma
    la t0, re_scratch
    li t1, 4
    /* Execute 0x1A628027 — should SIGILL */
    .4byte 0x1A628027
    /* If we reach here, instruction did NOT trap — exit child 0 */
    SYS_EXIT 0

re_smew_parent:
    mv s6, a0
    la s7, re_smew_wstatus
    SYS_WAIT4 s6, s7

    SET_TEST_NUM 2
    la t0, re_smew_wstatus
    lw t0, 0(t0)
    andi t0, t0, 0x7f
    li t1, 4
    FAIL_IF_NE t0, t1

    /* SIGILL check: unit-stride load with reserved lumop=00001 */
    /* Encoding 0x02128407 — must trap as illegal instruction */
    SET_TEST_NUM 3
    SYS_CLONE
    bltz a0, re_lumop_clone_fail
    beqz a0, re_lumop_child
    j re_lumop_parent
re_lumop_clone_fail:
    FAIL_TEST

re_lumop_child:
    li t0, 4
    vsetvli t0, t0, e8, m1, ta, ma
    la t0, re_scratch
    li t1, 4
    /* Execute 0x02128407 — should SIGILL */
    .4byte 0x02128407
    /* If we reach here, instruction did NOT trap — exit child 0 */
    SYS_EXIT 0

re_lumop_parent:
    mv s6, a0
    la s7, re_lumop_wstatus
    SYS_WAIT4 s6, s7

    SET_TEST_NUM 4
    la t0, re_lumop_wstatus
    lw t0, 0(t0)
    andi t0, t0, 0x7f
    li t1, 4
    FAIL_IF_NE t0, t1

    /* SIGILL check: mew=1 indexed load */
    /* Encoding 0x17028407 — must trap as illegal instruction */
    SET_TEST_NUM 5
    SYS_CLONE
    bltz a0, re_imew_clone_fail
    beqz a0, re_imew_child
    j re_imew_parent
re_imew_clone_fail:
    FAIL_TEST

re_imew_child:
    li t0, 4
    vsetvli t0, t0, e8, m1, ta, ma
    la t0, re_scratch
    li t1, 4
    /* Execute 0x17028407 — should SIGILL */
    .4byte 0x17028407
    /* If we reach here, instruction did NOT trap — exit child 0 */
    SYS_EXIT 0

re_imew_parent:
    mv s6, a0
    la s7, re_imew_wstatus
    SYS_WAIT4 s6, s7

    SET_TEST_NUM 6
    la t0, re_imew_wstatus
    lw t0, 0(t0)
    andi t0, t0, 0x7f
    li t1, 4
    FAIL_IF_NE t0, t1

    /* SIGILL check: unit-stride store with reserved sumop=00001 */
    /* Encoding 0x02128427 — must trap as illegal instruction */
    SET_TEST_NUM 7
    SYS_CLONE
    bltz a0, re_sumop_clone_fail
    beqz a0, re_sumop_child
    j re_sumop_parent
re_sumop_clone_fail:
    FAIL_TEST

re_sumop_child:
    li t0, 4
    vsetvli t0, t0, e8, m1, ta, ma
    la t0, re_scratch
    li t1, 4
    /* Execute 0x02128427 — should SIGILL */
    .4byte 0x02128427
    /* If we reach here, instruction did NOT trap — exit child 0 */
    SYS_EXIT 0

re_sumop_parent:
    mv s6, a0
    la s7, re_sumop_wstatus
    SYS_WAIT4 s6, s7

    SET_TEST_NUM 8
    la t0, re_sumop_wstatus
    lw t0, 0(t0)
    andi t0, t0, 0x7f
    li t1, 4
    FAIL_IF_NE t0, t1

    PASS_TEST

.data
.align 4
re_scratch:
    .space 64
.align 3
re_smew_wstatus:
    .word 0
re_lumop_wstatus:
    .word 0
re_imew_wstatus:
    .word 0
re_sumop_wstatus:
    .word 0

.align 4
result_buf:  .space 256
witness_buf: .space 256

