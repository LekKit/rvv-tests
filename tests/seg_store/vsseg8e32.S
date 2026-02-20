/* Auto-generated test for vsseg8e32.v
 * Unit-stride segment store nf=8 eew=32
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsseg8e32.v: memory result
 *     2 = vsseg8e32.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc1_f0
    vle32.v v8, (t1)
    la t1, tc1_f1
    vle32.v v9, (t1)
    la t1, tc1_f2
    vle32.v v10, (t1)
    la t1, tc1_f3
    vle32.v v11, (t1)
    la t1, tc1_f4
    vle32.v v12, (t1)
    la t1, tc1_f5
    vle32.v v13, (t1)
    la t1, tc1_f6
    vle32.v v14, (t1)
    la t1, tc1_f7
    vle32.v v15, (t1)
    SAVE_CSRS
    la t1, result_buf
    vsseg8e32.v v8, (t1)
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 128
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tc1_f0:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc1_f1:
    .word 0x00000005, 0x00000006, 0x00000007, 0x00000008
tc1_f2:
    .word 0x00000009, 0x0000000a, 0x0000000b, 0x0000000c
tc1_f3:
    .word 0x0000000d, 0x0000000e, 0x0000000f, 0x00000010
tc1_f4:
    .word 0x00000011, 0x00000012, 0x00000013, 0x00000014
tc1_f5:
    .word 0x00000015, 0x00000016, 0x00000017, 0x00000018
tc1_f6:
    .word 0x00000019, 0x0000001a, 0x0000001b, 0x0000001c
tc1_f7:
    .word 0x0000001d, 0x0000001e, 0x0000001f, 0x00000020
tc1_exp:
    .word 0x00000001, 0x00000005, 0x00000009, 0x0000000d, 0x00000011, 0x00000015, 0x00000019, 0x0000001d, 0x00000002, 0x00000006, 0x0000000a, 0x0000000e, 0x00000012, 0x00000016, 0x0000001a, 0x0000001e, 0x00000003, 0x00000007, 0x0000000b, 0x0000000f, 0x00000013, 0x00000017, 0x0000001b, 0x0000001f, 0x00000004, 0x00000008, 0x0000000c, 0x00000010, 0x00000014, 0x00000018, 0x0000001c, 0x00000020

.align 4
result_buf:  .space 256
witness_buf: .space 256

