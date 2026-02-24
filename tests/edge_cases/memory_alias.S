/* Auto-generated test for Memory aliasing coherence
 * Tests vector load/store coherence when two virtual addresses map to the same physical page (memfd + double mmap). Catches store buffer and cache coherence bugs.
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = memfd_create failed
 *     2 = ftruncate failed
 *     3 = mmap VA1 failed
 *     4 = mmap VA2 failed
 *     5 = vse through VA1 → vle from VA2: data mismatch
 *     6 = scalar sw through VA1 → vle from VA2: data mismatch
 *     7 = vse through VA1 → scalar lw from VA2: data mismatch
 *     8 = e8 alias: vse8 through VA1, vle8 from VA2 mismatch
 *     9 = e64 alias: vse64 through VA1, vle64 from VA2 mismatch
 *    10 = strided alias: vsse32 through VA1, vle32 from VA2 mismatch
 *    11 = alias loop: iteration mismatch
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Memory aliasing: two VAs → same physical page */
    /* Setup: memfd_create + ftruncate + mmap × 2 */
    SET_TEST_NUM 1
    la t0, malias_name
    SYS_MEMFD_CREATE t0
    bltz a0, malias_fail
    mv s6, a0
    SET_TEST_NUM 2
    SYS_FTRUNCATE s6, 4096
    bltz a0, malias_fail
    SET_TEST_NUM 3
    SYS_MMAP_FD 0, 4096, 3, 1, s6, 0
    bltz a0, malias_fail
    mv s7, a0
    SET_TEST_NUM 4
    SYS_MMAP_FD 0, 4096, 3, 1, s6, 0
    bltz a0, malias_fail
    mv s8, a0
    SYS_CLOSE s6

    /* Test 1: vector store through VA1, vector load from VA2 */
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    la t1, malias_data1
    vle32.v v16, (t1)
    vse32.v v16, (s7)
    fence rw, rw
    vle32.v v8, (s8)
    SET_TEST_NUM 5
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, malias_data1, 16

    /* Test 2: scalar stores through VA1, vector load from VA2 */
    la t1, malias_data2
    lw t2, 0(t1)
    sw t2, 0(s7)
    lw t2, 4(t1)
    sw t2, 4(s7)
    lw t2, 8(t1)
    sw t2, 8(s7)
    lw t2, 12(t1)
    sw t2, 12(s7)
    fence rw, rw
    vle32.v v8, (s8)
    SET_TEST_NUM 6
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, malias_data2, 16

    /* Test 3: vector store through VA1, scalar loads from VA2 */
    la t1, malias_data3
    vle32.v v16, (t1)
    vse32.v v16, (s7)
    fence rw, rw
    SET_TEST_NUM 7
    lw t2, 0(s8)
    la t1, malias_data3
    lw t3, 0(t1)
    FAIL_IF_NE t2, t3
    lw t2, 4(s8)
    la t1, malias_data3
    lw t3, 4(t1)
    FAIL_IF_NE t2, t3
    lw t2, 8(s8)
    la t1, malias_data3
    lw t3, 8(t1)
    FAIL_IF_NE t2, t3
    lw t2, 12(s8)
    la t1, malias_data3
    lw t3, 12(t1)
    FAIL_IF_NE t2, t3

    /* Test 4: e8 aliased store/load */
    li t0, 4
    vsetvli t0, t0, e8, m1, ta, ma
    la t1, malias_data_e8
    vle8.v v16, (t1)
    vse8.v v16, (s7)
    fence rw, rw
    vle8.v v8, (s8)
    SET_TEST_NUM 8
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, malias_data_e8, 4

    /* Test 5: e64 aliased store/load */
    li t0, 4
    vsetvli t0, t0, e64, m1, ta, ma
    la t1, malias_data_e64
    vle64.v v16, (t1)
    vse64.v v16, (s7)
    fence rw, rw
    vle64.v v8, (s8)
    SET_TEST_NUM 9
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, malias_data_e64, 32

    /* Test 6: strided store through alias, unit load back */
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    la t1, malias_data1
    vle32.v v16, (t1)
    li t2, 8
    vsse32.v v16, (s7), t2
    fence rw, rw
    vlse32.v v8, (s8), t2
    SET_TEST_NUM 10
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, malias_data1, 16

    /* Test 7: 4 write-read iterations through alias */
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    li s9, 0
malias_loop:
    slli t2, s9, 2
    vmv.v.x v16, t2
    vid.v v20
    vadd.vv v16, v16, v20
    vse32.v v16, (s7)
    fence rw, rw
    vle32.v v8, (s8)
    vmseq.vv v0, v8, v16
    vcpop.m t3, v0
    li t4, 4
    SET_TEST_NUM 11
    bne t3, t4, malias_fail
    addi s9, s9, 1
    li t4, 4
    blt s9, t4, malias_loop

    /* Cleanup */
    SYS_MUNMAP s7, 4096
    SYS_MUNMAP s8, 4096
    j malias_done
malias_fail:
    FAIL_TEST
malias_done:

    PASS_TEST

.data
.align 2
malias_data1:
    .word 0x11111111, 0x22222222, 0x33333333, 0x44444444
malias_data2:
    .word 0xaaaabbbb, 0xccccdddd, 0xeeeeffff, 0x12345678
malias_data3:
    .word 0xdeadbeef, 0xcafebabe, 0xfacefeed, 0xbaadf00d
.align 1
malias_data_e8:
    .byte 0xaa, 0xbb, 0xcc, 0xdd
.align 3
malias_data_e64:
    .dword 0x1111111122222222, 0x3333333344444444, 0x5555555566666666, 0x7777777788888888
malias_name:
    .asciz "rvvtest"

.align 4
result_buf:  .space 256
witness_buf: .space 256

