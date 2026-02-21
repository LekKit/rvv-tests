/* Auto-generated test for vfadd.vv
 * LMUL>1 register boundary crossing: vfadd.vv e32/m2 — 16 elements spanning v8-v9 on VLEN=256
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = lmul>1 fp: wrong result
 *     2 = lmul>1 fp: witness changed
 *     3 = lmul>1 fp: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* LMUL>1 test: vfadd.vv e32/m2, vl=16 */
    li t0, 16
    vsetvli t0, t0, e32, m2, tu, mu
    la t1, tc1_s2
    vle32.v v16, (t1)
    la t1, tc1_s1
    vle32.v v20, (t1)
    la t1, tc1_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vfadd.vv v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 64
    SET_TEST_NUM 2
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 64
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED_FP

    PASS_TEST

.data
.align 2
tc1_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000, 0x40a00000, 0x40c00000, 0x40e00000, 0x41000000, 0x41100000, 0x41200000, 0x41300000, 0x41400000, 0x41500000, 0x41600000, 0x41700000, 0x41800000
tc1_s1:
    .word 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000
tc1_exp:
    .word 0x3fc00000, 0x40200000, 0x40600000, 0x40900000, 0x40b00000, 0x40d00000, 0x40f00000, 0x41080000, 0x41180000, 0x41280000, 0x41380000, 0x41480000, 0x41580000, 0x41680000, 0x41780000, 0x41840000
tc1_w:
    .word 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef

.align 4
result_buf:  .space 256
witness_buf: .space 256

