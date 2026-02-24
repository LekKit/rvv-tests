/* Auto-generated test for Whole register operations
 * Tests vmv2r/4r/8r move data integrity and vl2re/vl4re/vs2r/vs4r load/store round-trips
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vmv2r.v: data integrity
 *     2 = vmv4r.v: data integrity
 *     3 = vl2re32/vs2r round-trip: data integrity
 *     4 = vl4re32/vs4r round-trip: data integrity
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* vmv2r.v: copy 2-register group, verify via load+move+store */
    la t1, wreg1_mem
    vl2re32.v v16, (t1)
    vmv2r.v v8, v16
    la t1, wreg1_res
    vs2r.v v8, (t1)
    SET_TEST_NUM 1
    csrr a3, vlenb
    slli a3, a3, 1
    la a1, wreg1_res
    la a2, wreg1_mem
    call _mem_compare
    FAIL_IF_NZ a0

    /* vmv4r.v: copy 4-register group */
    la t1, wreg2_mem
    vl4re32.v v16, (t1)
    vmv4r.v v8, v16
    la t1, wreg2_res
    vs4r.v v8, (t1)
    SET_TEST_NUM 2
    csrr a3, vlenb
    slli a3, a3, 2
    la a1, wreg2_res
    la a2, wreg2_mem
    call _mem_compare
    FAIL_IF_NZ a0

    /* vl2re32.v + vs2r.v round-trip */
    la t1, wreg3_mem
    vl2re32.v v8, (t1)
    la t1, wreg3_res
    vs2r.v v8, (t1)
    SET_TEST_NUM 3
    csrr a3, vlenb
    slli a3, a3, 1
    la a1, wreg3_res
    la a2, wreg3_mem
    call _mem_compare
    FAIL_IF_NZ a0

    /* vl4re32.v + vs4r.v round-trip */
    la t1, wreg4_mem
    vl4re32.v v8, (t1)
    la t1, wreg4_res
    vs4r.v v8, (t1)
    SET_TEST_NUM 4
    csrr a3, vlenb
    slli a3, a3, 2
    la a1, wreg4_res
    la a2, wreg4_mem
    call _mem_compare
    FAIL_IF_NZ a0

    PASS_TEST

.data
.align 4
wreg1_mem:
    .word 0x11111111, 0x12121212, 0x13131313, 0x14141414, 0x15151515, 0x16161616, 0x17171717, 0x18181818, 0x19191919, 0x1a1a1a1a, 0x1b1b1b1b, 0x1c1c1c1c, 0x1d1d1d1d, 0x1e1e1e1e, 0x1f1f1f1f, 0x20202020, 0x21212121, 0x22222222, 0x23232323, 0x24242424, 0x25252525, 0x26262626, 0x27272727, 0x28282828, 0x29292929, 0x2a2a2a2a, 0x2b2b2b2b, 0x2c2c2c2c, 0x2d2d2d2d, 0x2e2e2e2e, 0x2f2f2f2f, 0x30303030
wreg1_res:
    .space 128
.align 4
wreg2_mem:
    .word 0xaaaa0001, 0xaaaa0002, 0xaaaa0003, 0xaaaa0004, 0xaaaa0005, 0xaaaa0006, 0xaaaa0007, 0xaaaa0008, 0xaaaa0009, 0xaaaa000a, 0xaaaa000b, 0xaaaa000c, 0xaaaa000d, 0xaaaa000e, 0xaaaa000f, 0xaaaa0010, 0xaaaa0011, 0xaaaa0012, 0xaaaa0013, 0xaaaa0014, 0xaaaa0015, 0xaaaa0016, 0xaaaa0017, 0xaaaa0018, 0xaaaa0019, 0xaaaa001a, 0xaaaa001b, 0xaaaa001c, 0xaaaa001d, 0xaaaa001e, 0xaaaa001f, 0xaaaa0020, 0xaaaa0021, 0xaaaa0022, 0xaaaa0023, 0xaaaa0024, 0xaaaa0025, 0xaaaa0026, 0xaaaa0027, 0xaaaa0028, 0xaaaa0029, 0xaaaa002a, 0xaaaa002b, 0xaaaa002c, 0xaaaa002d, 0xaaaa002e, 0xaaaa002f, 0xaaaa0030, 0xaaaa0031, 0xaaaa0032, 0xaaaa0033, 0xaaaa0034, 0xaaaa0035, 0xaaaa0036, 0xaaaa0037, 0xaaaa0038, 0xaaaa0039, 0xaaaa003a, 0xaaaa003b, 0xaaaa003c, 0xaaaa003d, 0xaaaa003e, 0xaaaa003f, 0xaaaa0040
wreg2_res:
    .space 256
.align 4
wreg3_mem:
    .word 0x12345678, 0x9abcdef0, 0xfedcba98, 0x76543210
    .space 112
wreg3_res:
    .space 128
.align 4
wreg4_mem:
    .word 0x12345678, 0x9abcdef0, 0xfedcba98, 0x76543210
    .space 240
wreg4_res:
    .space 256

.align 4
result_buf:  .space 256
witness_buf: .space 256

