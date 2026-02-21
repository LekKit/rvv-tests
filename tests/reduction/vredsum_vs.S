/* Auto-generated test for vredsum.vs
 * Reduction vredsum.vs
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vredsum.vs e8 basic: wrong result
 *     2 = vredsum.vs e8 basic: witness changed
 *     3 = vredsum.vs e8 basic: CSR side-effect
 *     4 = vredsum.vs e8 init: wrong result
 *     5 = vredsum.vs e8 init: witness changed
 *     6 = vredsum.vs e8 init: CSR side-effect
 *     7 = vredsum.vs e8 max: wrong result
 *     8 = vredsum.vs e8 max: witness changed
 *     9 = vredsum.vs e8 max: CSR side-effect
 *    10 = vredsum.vs e16 basic: wrong result
 *    11 = vredsum.vs e16 basic: witness changed
 *    12 = vredsum.vs e16 basic: CSR side-effect
 *    13 = vredsum.vs e16 init: wrong result
 *    14 = vredsum.vs e16 init: witness changed
 *    15 = vredsum.vs e16 init: CSR side-effect
 *    16 = vredsum.vs e16 max: wrong result
 *    17 = vredsum.vs e16 max: witness changed
 *    18 = vredsum.vs e16 max: CSR side-effect
 *    19 = vredsum.vs e32 basic: wrong result
 *    20 = vredsum.vs e32 basic: witness changed
 *    21 = vredsum.vs e32 basic: CSR side-effect
 *    22 = vredsum.vs e32 init: wrong result
 *    23 = vredsum.vs e32 init: witness changed
 *    24 = vredsum.vs e32 init: CSR side-effect
 *    25 = vredsum.vs e32 max: wrong result
 *    26 = vredsum.vs e32 max: witness changed
 *    27 = vredsum.vs e32 max: CSR side-effect
 *    28 = vredsum.vs e64 basic: wrong result
 *    29 = vredsum.vs e64 basic: witness changed
 *    30 = vredsum.vs e64 basic: CSR side-effect
 *    31 = vredsum.vs e64 init: wrong result
 *    32 = vredsum.vs e64 init: witness changed
 *    33 = vredsum.vs e64 init: CSR side-effect
 *    34 = vredsum.vs e64 max: wrong result
 *    35 = vredsum.vs e64 max: witness changed
 *    36 = vredsum.vs e64 max: CSR side-effect
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Test 1-3: vredsum.vs SEW=8 basic */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc1_vec
    vle8.v v16, (t1)
    la t1, tc1_init
    vle8.v v20, (t1)
    la t1, tc1_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 1
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc1_exp, 1
    SET_TEST_NUM 2
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc1_w, 4
    SET_TEST_NUM 3
    CHECK_CSRS_UNCHANGED

    /* Test 4-6: vredsum.vs SEW=8 init */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc4_vec
    vle8.v v16, (t1)
    la t1, tc4_init
    vle8.v v20, (t1)
    la t1, tc4_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 4
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc4_exp, 1
    SET_TEST_NUM 5
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc4_w, 4
    SET_TEST_NUM 6
    CHECK_CSRS_UNCHANGED

    /* Test 7-9: vredsum.vs SEW=8 max */
    li t0, 4
    vsetvli t0, t0, e8, m1, tu, mu
    la t1, tc7_vec
    vle8.v v16, (t1)
    la t1, tc7_init
    vle8.v v20, (t1)
    la t1, tc7_w
    vle8.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 7
    la t1, result_buf
    vse8.v v8, (t1)
    CHECK_MEM result_buf, tc7_exp, 1
    SET_TEST_NUM 8
    la t1, witness_buf
    vse8.v v24, (t1)
    CHECK_MEM witness_buf, tc7_w, 4
    SET_TEST_NUM 9
    CHECK_CSRS_UNCHANGED

    /* Test 10-12: vredsum.vs SEW=16 basic */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc10_vec
    vle16.v v16, (t1)
    la t1, tc10_init
    vle16.v v20, (t1)
    la t1, tc10_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 10
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc10_exp, 2
    SET_TEST_NUM 11
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc10_w, 8
    SET_TEST_NUM 12
    CHECK_CSRS_UNCHANGED

    /* Test 13-15: vredsum.vs SEW=16 init */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc13_vec
    vle16.v v16, (t1)
    la t1, tc13_init
    vle16.v v20, (t1)
    la t1, tc13_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 13
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc13_exp, 2
    SET_TEST_NUM 14
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc13_w, 8
    SET_TEST_NUM 15
    CHECK_CSRS_UNCHANGED

    /* Test 16-18: vredsum.vs SEW=16 max */
    li t0, 4
    vsetvli t0, t0, e16, m1, tu, mu
    la t1, tc16_vec
    vle16.v v16, (t1)
    la t1, tc16_init
    vle16.v v20, (t1)
    la t1, tc16_w
    vle16.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 16
    la t1, result_buf
    vse16.v v8, (t1)
    CHECK_MEM result_buf, tc16_exp, 2
    SET_TEST_NUM 17
    la t1, witness_buf
    vse16.v v24, (t1)
    CHECK_MEM witness_buf, tc16_w, 8
    SET_TEST_NUM 18
    CHECK_CSRS_UNCHANGED

    /* Test 19-21: vredsum.vs SEW=32 basic */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc19_vec
    vle32.v v16, (t1)
    la t1, tc19_init
    vle32.v v20, (t1)
    la t1, tc19_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 19
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc19_exp, 4
    SET_TEST_NUM 20
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc19_w, 16
    SET_TEST_NUM 21
    CHECK_CSRS_UNCHANGED

    /* Test 22-24: vredsum.vs SEW=32 init */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc22_vec
    vle32.v v16, (t1)
    la t1, tc22_init
    vle32.v v20, (t1)
    la t1, tc22_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 22
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc22_exp, 4
    SET_TEST_NUM 23
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc22_w, 16
    SET_TEST_NUM 24
    CHECK_CSRS_UNCHANGED

    /* Test 25-27: vredsum.vs SEW=32 max */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tc25_vec
    vle32.v v16, (t1)
    la t1, tc25_init
    vle32.v v20, (t1)
    la t1, tc25_w
    vle32.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 25
    la t1, result_buf
    vse32.v v8, (t1)
    CHECK_MEM result_buf, tc25_exp, 4
    SET_TEST_NUM 26
    la t1, witness_buf
    vse32.v v24, (t1)
    CHECK_MEM witness_buf, tc25_w, 16
    SET_TEST_NUM 27
    CHECK_CSRS_UNCHANGED

    /* Test 28-30: vredsum.vs SEW=64 basic */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc28_vec
    vle64.v v16, (t1)
    la t1, tc28_init
    vle64.v v20, (t1)
    la t1, tc28_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 28
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc28_exp, 8
    SET_TEST_NUM 29
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc28_w, 32
    SET_TEST_NUM 30
    CHECK_CSRS_UNCHANGED

    /* Test 31-33: vredsum.vs SEW=64 init */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc31_vec
    vle64.v v16, (t1)
    la t1, tc31_init
    vle64.v v20, (t1)
    la t1, tc31_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 31
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc31_exp, 8
    SET_TEST_NUM 32
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc31_w, 32
    SET_TEST_NUM 33
    CHECK_CSRS_UNCHANGED

    /* Test 34-36: vredsum.vs SEW=64 max */
    li t0, 4
    vsetvli t0, t0, e64, m1, tu, mu
    la t1, tc34_vec
    vle64.v v16, (t1)
    la t1, tc34_init
    vle64.v v20, (t1)
    la t1, tc34_w
    vle64.v v24, (t1)
    SAVE_CSRS
    vredsum.vs v8, v16, v20
    SET_TEST_NUM 34
    la t1, result_buf
    vse64.v v8, (t1)
    CHECK_MEM result_buf, tc34_exp, 8
    SET_TEST_NUM 35
    la t1, witness_buf
    vse64.v v24, (t1)
    CHECK_MEM witness_buf, tc34_w, 32
    SET_TEST_NUM 36
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 1
tc1_vec:
    .byte 0x01, 0x02, 0x03, 0x04
tc1_init:
    .byte 0x00, 0x00, 0x00, 0x00
tc1_exp:
    .byte 0x0a, 0x00, 0x00, 0x00
tc1_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc4_vec:
    .byte 0x01, 0x02, 0x03, 0x04
tc4_init:
    .byte 0x64, 0x00, 0x00, 0x00
tc4_exp:
    .byte 0x6e, 0x00, 0x00, 0x00
tc4_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc7_vec:
    .byte 0xff, 0xfe, 0xfd, 0x00
tc7_init:
    .byte 0x00, 0x00, 0x00, 0x00
tc7_exp:
    .byte 0xfa, 0x00, 0x00, 0x00
tc7_w:
    .byte 0xde, 0xad, 0xbe, 0xef
.align 1
tc10_vec:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc10_init:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc10_exp:
    .half 0x000a, 0x0000, 0x0000, 0x0000
tc10_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc13_vec:
    .half 0x0001, 0x0002, 0x0003, 0x0004
tc13_init:
    .half 0x0064, 0x0000, 0x0000, 0x0000
tc13_exp:
    .half 0x006e, 0x0000, 0x0000, 0x0000
tc13_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 1
tc16_vec:
    .half 0xffff, 0xfffe, 0xfffd, 0x0000
tc16_init:
    .half 0x0000, 0x0000, 0x0000, 0x0000
tc16_exp:
    .half 0xfffa, 0x0000, 0x0000, 0x0000
tc16_w:
    .half 0xdead, 0xbeef, 0xcafe, 0xbabe
.align 2
tc19_vec:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc19_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc19_exp:
    .word 0x0000000a, 0x00000000, 0x00000000, 0x00000000
tc19_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc22_vec:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tc22_init:
    .word 0x00000064, 0x00000000, 0x00000000, 0x00000000
tc22_exp:
    .word 0x0000006e, 0x00000000, 0x00000000, 0x00000000
tc22_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 2
tc25_vec:
    .word 0xffffffff, 0xfffffffe, 0xfffffffd, 0x00000000
tc25_init:
    .word 0x00000000, 0x00000000, 0x00000000, 0x00000000
tc25_exp:
    .word 0xfffffffa, 0x00000000, 0x00000000, 0x00000000
tc25_w:
    .word 0xdeadbeef, 0xcafebabe, 0x12345678, 0x9abcdef0
.align 3
tc28_vec:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc28_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc28_exp:
    .dword 0x000000000000000a, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc28_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc31_vec:
    .dword 0x0000000000000001, 0x0000000000000002, 0x0000000000000003, 0x0000000000000004
tc31_init:
    .dword 0x0000000000000064, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc31_exp:
    .dword 0x000000000000006e, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc31_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef
.align 3
tc34_vec:
    .dword 0xffffffffffffffff, 0xfffffffffffffffe, 0xfffffffffffffffd, 0x0000000000000000
tc34_init:
    .dword 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc34_exp:
    .dword 0xfffffffffffffffa, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000
tc34_w:
    .dword 0xdeadbeefcafebabe, 0x123456789abcdef0, 0xfedcba9876543210, 0x0123456789abcdef

.align 4
result_buf:  .space 256
witness_buf: .space 256

