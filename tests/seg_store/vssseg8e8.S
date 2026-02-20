/* Auto-generated test for vssseg8e8.v
 * Strided segment store nf=8 eew=8
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vssseg8e8.v: memory result
 *     2 = vssseg8e8.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_f0
    vle8.v v8, (t1)
    la t1, tc1_f1
    vle8.v v9, (t1)
    la t1, tc1_f2
    vle8.v v10, (t1)
    la t1, tc1_f3
    vle8.v v11, (t1)
    la t1, tc1_f4
    vle8.v v12, (t1)
    la t1, tc1_f5
    vle8.v v13, (t1)
    la t1, tc1_f6
    vle8.v v14, (t1)
    la t1, tc1_f7
    vle8.v v15, (t1)
    SAVE_CSRS
    la t1, result_buf
    li t2, 8
    vssseg8e8.v v8, (t1), t2
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 32
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_f0:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_f1:
    .byte 0x05, 0x06, 0x07, 0x08
tc1_f2:
    .byte 0x09, 0x0a, 0x0b, 0x0c
tc1_f3:
    .byte 0x0d, 0x0e, 0x0f, 0x10
tc1_f4:
    .byte 0x11, 0x12, 0x13, 0x14
tc1_f5:
    .byte 0x15, 0x16, 0x17, 0x18
tc1_f6:
    .byte 0x19, 0x1a, 0x1b, 0x1c
tc1_f7:
    .byte 0x1d, 0x1e, 0x1f, 0x20
tc1_exp:
    .byte 0x01, 0x05, 0x09, 0x0d, 0x11, 0x15, 0x19, 0x1d, 0x02, 0x06, 0x0a, 0x0e, 0x12, 0x16, 0x1a, 0x1e, 0x03, 0x07, 0x0b, 0x0f, 0x13, 0x17, 0x1b, 0x1f, 0x04, 0x08, 0x0c, 0x10, 0x14, 0x18, 0x1c, 0x20

.align 4
result_buf:  .space 256
witness_buf: .space 256

