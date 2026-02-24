/* Auto-generated test for vloxei64.v
 * Indexed load vloxei64.v
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vloxei64.v identity: result
 *     2 = vloxei64.v reverse: result
 *     3 = vloxei64.v splat: result
 *     4 = vloxei64.v masked: active+inactive
 */
#include "riscv_test.h"
#include "test_macros.h"

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc1_idx
    vle64.v v20, (t1)
    la t1, tc1_data
    SAVE_CSRS
    vloxei64.v v8, (t1), v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc2_idx
    vle64.v v20, (t1)
    la t1, tc2_data
    SAVE_CSRS
    vloxei64.v v8, (t1), v20
    SET_TEST_NUM 2
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc2_exp, 32
    CHECK_CSRS_UNCHANGED

    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc3_idx
    vle64.v v20, (t1)
    la t1, tc3_data
    SAVE_CSRS
    vloxei64.v v8, (t1), v20
    SET_TEST_NUM 3
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc3_exp, 32
    CHECK_CSRS_UNCHANGED

    /* Masked indexed load: mask=0b0101 */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc4_vd
    vle64.v v8, (t1)
    la t1, tc4_idx
    vle64.v v20, (t1)
    la t1, tc4_mask
    vlm.v v0, (t1)
    SAVE_CSRS
    la t1, tc4_data
    vloxei64.v v8, (t1), v20, v0.t
    SET_TEST_NUM 4
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 32
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 3
tc1_data:
    .dword 0x000000000000000a, 0x000000000000000b, 0x000000000000000c, 0x000000000000000d, 0x000000000000000e, 0x000000000000000f, 0x0000000000000010, 0x0000000000000011
tc1_idx:
    .dword 0x0000000000000000, 0x0000000000000008, 0x0000000000000010, 0x0000000000000018
tc1_exp:
    .dword 0x000000000000000a, 0x000000000000000b, 0x000000000000000c, 0x000000000000000d
.align 3
tc2_data:
    .dword 0x000000000000000a, 0x000000000000000b, 0x000000000000000c, 0x000000000000000d, 0x000000000000000e, 0x000000000000000f, 0x0000000000000010, 0x0000000000000011
tc2_idx:
    .dword 0x0000000000000018, 0x0000000000000010, 0x0000000000000008, 0x0000000000000000
tc2_exp:
    .dword 0x000000000000000d, 0x000000000000000c, 0x000000000000000b, 0x000000000000000a
.align 3
tc3_data:
    .dword 0x000000000000000a, 0x000000000000000b, 0x000000000000000c, 0x000000000000000d, 0x000000000000000e, 0x000000000000000f, 0x0000000000000010, 0x0000000000000011
tc3_idx:
    .dword 0x0000000000000010, 0x0000000000000010, 0x0000000000000010, 0x0000000000000010
tc3_exp:
    .dword 0x000000000000000c, 0x000000000000000c, 0x000000000000000c, 0x000000000000000c
.align 1
tc4_mask:
    .byte 5
.align 3
tc4_vd:
    .dword 0x00000000000000dd, 0x00000000000000dd, 0x00000000000000dd, 0x00000000000000dd
tc4_data:
    .dword 0x000000000000000a, 0x000000000000000b, 0x000000000000000c, 0x000000000000000d, 0x000000000000000e, 0x000000000000000f, 0x0000000000000010, 0x0000000000000011
tc4_idx:
    .dword 0x0000000000000000, 0x0000000000000008, 0x0000000000000010, 0x0000000000000018
tc4_exp:
    .dword 0x000000000000000a, 0x00000000000000dd, 0x000000000000000c, 0x00000000000000dd

.align 4
result_buf:  .space 256
witness_buf: .space 256

