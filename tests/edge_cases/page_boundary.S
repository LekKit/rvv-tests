/* Auto-generated test for Page boundary vector access
 * Tests vector load/store spanning a page boundary (elements on page 1 and page 2). Verifies TLB refill mid-vector-op.
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = mmap 2 pages failed
 *     2 = cross-page vector store failed
 *     3 = cross-page vector load mismatch
 *     4 = cross-page scalar verify mismatch
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Page boundary: vector access spanning two pages */
    SET_TEST_NUM 1
    SYS_MMAP 0, 8192, 3, 34, -1, 0
    bltz a0, pgb_fail
    mv s6, a0
    li t0, 4088
    add s7, s6, t0
    SET_TEST_NUM 2
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    la t1, pgb_data
    vle32.v v16, (t1)
    vse32.v v16, (s7)
    vle32.v v8, (s7)
    SET_TEST_NUM 3
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, pgb_data, 16
    SET_TEST_NUM 4
    lw t2, 0(s7)
    la t1, pgb_data
    lw t3, 0(t1)
    FAIL_IF_NE t2, t3
    lw t2, 8(s7)
    lw t3, 8(t1)
    FAIL_IF_NE t2, t3
    SYS_MUNMAP s6, 8192
    j pgb_done
pgb_fail:
    FAIL_TEST
pgb_done:

    PASS_TEST

.data
.align 4
pgb_data:
    .word 0xaaaa0001, 0xbbbb0002, 0xcccc0003, 0xdddd0004

.align 4
result_buf:  .space 256
witness_buf: .space 256

