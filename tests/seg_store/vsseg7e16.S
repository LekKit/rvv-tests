/* Auto-generated test for vsseg7e16.v
 * Unit-stride segment store nf=7 eew=16
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vsseg7e16.v: memory result
 *     2 = vsseg7e16.v: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc1_f0
    vle16.v v8, (t1)
    la t1, tc1_f1
    vle16.v v9, (t1)
    la t1, tc1_f2
    vle16.v v10, (t1)
    la t1, tc1_f3
    vle16.v v11, (t1)
    la t1, tc1_f4
    vle16.v v12, (t1)
    la t1, tc1_f5
    vle16.v v13, (t1)
    la t1, tc1_f6
    vle16.v v14, (t1)
    SAVE_CSRS
    la t1, result_buf
    vsseg7e16.v v8, (t1)
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tc1_exp, 56
    SET_TEST_NUM 2
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_f0:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc1_f1:
    .half 0x0005, 0x0006, 0x0007, 0x0008
tc1_f2:
    .half 0x0009, 0x000a, 0x000b, 0x000c
tc1_f3:
    .half 0x000d, 0x000e, 0x000f, 0x0010
tc1_f4:
    .half 0x0011, 0x0012, 0x0013, 0x0014
tc1_f5:
    .half 0x0015, 0x0016, 0x0017, 0x0018
tc1_f6:
    .half 0x0019, 0x001a, 0x001b, 0x001c
tc1_exp:
    .half 0x0001, 0x0005, 0x0009, 0x000d, 0x0011, 0x0015, 0x0019, 0x0002, 0x0006, 0x000a, 0x000e, 0x0012, 0x0016, 0x001a, 0x0003, 0x0007, 0x000b, 0x000f, 0x0013, 0x0017, 0x001b, 0x0004, 0x0008, 0x000c, 0x0010, 0x0014, 0x0018, 0x001c

.align 4
result_buf:  .space 256
witness_buf: .space 256

