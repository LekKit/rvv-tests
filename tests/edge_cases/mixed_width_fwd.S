/* Auto-generated test for Mixed-width store-load forwarding
 * Tests vector store at one SEW followed by vector load at a different SEW from the same address. Verifies store buffer handles cross-width forwarding correctly.
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vse32 then vle8: byte mismatch
 *     2 = vse8 then vle32: word mismatch
 *     3 = vse64 then vle16: halfword mismatch
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vse32 then vle8 from same address */
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    la t1, mwf_data32
    vle32.v v16, (t1)
    la s6, mwf_buf1
    vse32.v v16, (s6)
    li t0, 16
    vsetvli t0, t0, e8, m1, ta, ma
    vle8.v v8, (s6)
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, mwf_exp_e8, 16

    /* vse8 then vle32 from same address */
    li t0, 16
    vsetvli t0, t0, e8, m1, ta, ma
    la t1, mwf_data_bytes
    vle8.v v16, (t1)
    la s6, mwf_buf2
    vse8.v v16, (s6)
    li t0, 4
    vsetvli t0, t0, e32, m1, ta, ma
    vle32.v v8, (s6)
    SET_TEST_NUM 2
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, mwf_exp_e32, 16

    /* vse64 then vle16 from same address */
    li t0, 4
    vsetvli t0, t0, e64, m1, ta, ma
    la t1, mwf_data64
    vle64.v v16, (t1)
    la s6, mwf_buf3
    vse64.v v16, (s6)
    li t0, 16
    vsetvli t0, t0, e16, m1, ta, ma
    vle16.v v8, (s6)
    SET_TEST_NUM 3
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, mwf_exp_e16, 32

    PASS_TEST

.data
.align 4
mwf_data32:
    .word 0x04030201, 0x08070605, 0x0c0b0a09, 0x100f0e0d
mwf_exp_e8:
    .byte 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10
mwf_data_bytes:
    .byte 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff, 0x00
mwf_exp_e32:
    .word 0x44332211, 0x88776655, 0xccbbaa99, 0x00ffeedd
.align 3
mwf_data64:
    .dword 0x0004000300020001, 0x0008000700060005, 0x000c000b000a0009, 0x0010000f000e000d
.align 1
mwf_exp_e16:
    .half 0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0006, 0x0007, 0x0008, 0x0009, 0x000a, 0x000b, 0x000c, 0x000d, 0x000e, 0x000f, 0x0010
.align 4
mwf_buf1:
    .space 32
mwf_buf2:
    .space 32
mwf_buf3:
    .space 64

.align 4
result_buf:  .space 256
witness_buf: .space 256

