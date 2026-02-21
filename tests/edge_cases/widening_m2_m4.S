/* Auto-generated test for vwadd.vv
 * Widening m2→m4: vwadd.vv e32/m2 → e64/m4 spanning 4 registers on VLEN=256
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = widening m2→m4: wrong result
 *     2 = widening m2→m4: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Widening test: vwadd.vv e32/m2 → e64/m4, vl=16 */
    /* Dest v8-v11 (m4), src v16-v17 (m2) + v20-v21 (m2) */
    li t0, 16
    vsetvli s6, t0, e32, m2, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vwadd.vv v8, v16, v20
    mv t0, s6
    vsetvli t0, t0, e64, m4, tu, mu
    la t1, tc1_resbuf
    vse64.v v8, (t1)
    SET_TEST_NUM 1
    la a1, tc1_resbuf
    la a2, tc1_exp
    slli a3, s6, 3
    jal ra, _mem_compare
    FAIL_IF_NZ a0
    SET_TEST_NUM 2
    CHECK_VSTART_ZERO
    CHECK_VXSAT_UNCHANGED
    CHECK_VXRM_UNCHANGED
    CHECK_FCSR_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_s2:
    .word 0x00000001, 0x0000000b, 0x00000015, 0x0000001f, 0x00000029, 0x00000033, 0x0000003d, 0x00000047, 0x00000051, 0x0000005b, 0x00000065, 0x0000006f, 0x00000079, 0x00000083, 0x0000008d, 0x00000097
tc1_s1:
    .word 0x00000002, 0x00000007, 0x0000000c, 0x00000011, 0x00000016, 0x0000001b, 0x00000020, 0x00000025, 0x0000002a, 0x0000002f, 0x00000034, 0x00000039, 0x0000003e, 0x00000043, 0x00000048, 0x0000004d
.align 3
tc1_exp:
    .dword 0x0000000000000003, 0x0000000000000012, 0x0000000000000021, 0x0000000000000030, 0x000000000000003f, 0x000000000000004e, 0x000000000000005d, 0x000000000000006c, 0x000000000000007b, 0x000000000000008a, 0x0000000000000099, 0x00000000000000a8, 0x00000000000000b7, 0x00000000000000c6, 0x00000000000000d5, 0x00000000000000e4
.align 4
tc1_resbuf:
    .space 256

.align 4
result_buf:  .space 256
witness_buf: .space 256

