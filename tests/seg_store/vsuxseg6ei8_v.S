/* Auto-generated test for vsuxseg6ei8.v
 * Indexed segment store nf=6 eew=8
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsuxseg6ei8.v: memory result
 *     2 = vsuxseg6ei8.v: CSR side-effect
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
    la t1, tc1_idx
    vle8.v v16, (t1)
    SAVE_CSRS
    la t1, result_buf
    vsuxseg6ei8.v v8, (t1), v16
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 24
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
tc1_idx:
    .byte 0x00, 0x06, 0x0c, 0x12
tc1_exp:
    .byte 0x01, 0x05, 0x09, 0x0d, 0x11, 0x15, 0x02, 0x06, 0x0a, 0x0e, 0x12, 0x16, 0x03, 0x07, 0x0b, 0x0f, 0x13, 0x17, 0x04, 0x08, 0x0c, 0x10, 0x14, 0x18

.align 4
result_buf:  .space 256
witness_buf: .space 256

