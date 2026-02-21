/* Auto-generated test for vadd.vv
 * LMUL>1 register boundary crossing: vadd.vv e32/m2 — 16 elements spanning v8-v9 on VLEN=256
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = lmul>1 int: wrong result
 *     2 = lmul>1 int: witness changed
 *     3 = lmul>1 int: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* LMUL>1 test: vadd.vv e32/m2, vl=16 */
    /* On VLEN=256: elements 0-7 in v8, elements 8-15 in v9 */
    li t0, 16
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    la t1, tc1_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vadd.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 64
    SET_TEST_NUM 2
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 64
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_s2:
    .word 0x00000001, 0x0000000b, 0x00000015, 0x0000001f, 0x00000029, 0x00000033, 0x0000003d, 0x00000047, 0x00000051, 0x0000005b, 0x00000065, 0x0000006f, 0x00000079, 0x00000083, 0x0000008d, 0x00000097
tc1_s1:
    .word 0x00000002, 0x00000007, 0x0000000c, 0x00000011, 0x00000016, 0x0000001b, 0x00000020, 0x00000025, 0x0000002a, 0x0000002f, 0x00000034, 0x00000039, 0x0000003e, 0x00000043, 0x00000048, 0x0000004d
tc1_exp:
    .word 0x00000003, 0x00000012, 0x00000021, 0x00000030, 0x0000003f, 0x0000004e, 0x0000005d, 0x0000006c, 0x0000007b, 0x0000008a, 0x00000099, 0x000000a8, 0x000000b7, 0x000000c6, 0x000000d5, 0x000000e4
tc1_w:
    .word 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef

.align 4
result_buf:  .space 256
witness_buf: .space 256

