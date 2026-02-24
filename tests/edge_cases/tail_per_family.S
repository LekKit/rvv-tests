/* Auto-generated test for Tail undisturbed per family
 * Verifies tail elements are preserved (tu policy) for representative instructions from each family: vadd, vmul, vwadd, vmacc, vsadd, vfadd, vle32, vslideup
 *
 * Exit code 0 = PASS
 * Exit code N = check N failed:
 *     1 = vadd.vv tu: tail not preserved
 *     2 = vmul.vv tu: tail not preserved
 *     3 = vsadd.vv tu: tail not preserved
 *     4 = vfadd.vv tu: tail not preserved
 *     5 = vle32.v tu: tail not preserved
 *     6 = vslideup.vi tu: tail not preserved
 */
#include "riscv_test.h"
#include "test_macros.h"


    /* Tail test: vadd.vv (vl=2, tu policy) */
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tpf1_pf
    vle32.v v8, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tpf1_s2
    vle32.v v16, (t1)
    la t1, tpf1_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vadd.vv v8, v16, v20
    la t1, result_buf
    vs1r.v v8, (t1)
    SET_TEST_NUM 1
    CHECK_MEM result_buf, tpf1_exp, 16
    CHECK_CSRS_UNCHANGED

    /* Tail test: vmul.vv (vl=2, tu policy) */
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tpf2_pf
    vle32.v v8, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tpf2_s2
    vle32.v v16, (t1)
    la t1, tpf2_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vmul.vv v8, v16, v20
    la t1, result_buf
    vs1r.v v8, (t1)
    SET_TEST_NUM 2
    CHECK_MEM result_buf, tpf2_exp, 16
    CHECK_CSRS_UNCHANGED

    /* Tail test: vsadd.vv (vl=2, tu policy) */
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tpf3_pf
    vle32.v v8, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tpf3_s2
    vle32.v v16, (t1)
    la t1, tpf3_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vsadd.vv v8, v16, v20
    la t1, result_buf
    vs1r.v v8, (t1)
    SET_TEST_NUM 3
    CHECK_MEM result_buf, tpf3_exp, 16
    CHECK_CSRS_UNCHANGED_FIXEDPOINT

    /* Tail test: vfadd.vv (vl=2, tu policy) */
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tpf4_pf
    vle32.v v8, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tpf4_s2
    vle32.v v16, (t1)
    la t1, tpf4_s1
    vle32.v v20, (t1)
    SAVE_CSRS
    vfadd.vv v8, v16, v20
    la t1, result_buf
    vs1r.v v8, (t1)
    SET_TEST_NUM 4
    CHECK_MEM result_buf, tpf4_exp, 16
    CHECK_CSRS_UNCHANGED_FP

    /* Tail test: vle32.v (vl=2, tu policy) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tpf5_pf
    vle32.v v8, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    la t1, tpf5_mem
    vle32.v v8, (t1)
    la t1, result_buf
    vs1r.v v8, (t1)
    SET_TEST_NUM 5
    CHECK_MEM result_buf, tpf5_exp, 16
    CHECK_CSRS_UNCHANGED

    /* Tail test: vslideup.vi (vl=2, offset=1, tu policy) */
    li t0, 4
    vsetvli t0, t0, e32, m1, tu, mu
    la t1, tpf6_pf
    vle32.v v8, (t1)
    la t1, tpf6_src
    vle32.v v16, (t1)
    li t0, 2
    vsetvli t0, t0, e32, m1, tu, mu
    SAVE_CSRS
    vslideup.vi v8, v16, 1
    la t1, result_buf
    vs1r.v v8, (t1)
    SET_TEST_NUM 6
    CHECK_MEM result_buf, tpf6_exp, 16
    CHECK_CSRS_UNCHANGED

    PASS_TEST

.data
.align 2
tpf1_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
tpf1_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tpf1_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tpf1_exp:
    .word 0x0000000b, 0x00000016, 0xcccccccc, 0xdddddddd
.align 2
tpf2_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
tpf2_s2:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tpf2_s1:
    .word 0x00000001, 0x00000002, 0x00000003, 0x00000004
tpf2_exp:
    .word 0x0000000a, 0x00000028, 0xcccccccc, 0xdddddddd
.align 2
tpf3_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
tpf3_s2:
    .word 0x7ffffffe, 0x7fffffff, 0x0000000a, 0x00000014
tpf3_s1:
    .word 0x00000001, 0x00000001, 0x00000005, 0x00000005
tpf3_exp:
    .word 0x7fffffff, 0x7fffffff, 0xcccccccc, 0xdddddddd
.align 2
tpf4_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
tpf4_s2:
    .word 0x3f800000, 0x40000000, 0x40400000, 0x40800000
tpf4_s1:
    .word 0x3f000000, 0x3f000000, 0x3f000000, 0x3f000000
tpf4_exp:
    .word 0x3fc00000, 0x40200000, 0xcccccccc, 0xdddddddd
.align 2
tpf5_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
tpf5_mem:
    .word 0x00001111, 0x00002222, 0x00003333, 0x00004444
tpf5_exp:
    .word 0x00001111, 0x00002222, 0xcccccccc, 0xdddddddd
.align 2
tpf6_pf:
    .word 0xaaaaaaaa, 0xbbbbbbbb, 0xcccccccc, 0xdddddddd
tpf6_src:
    .word 0x0000000a, 0x00000014, 0x0000001e, 0x00000028
tpf6_exp:
    .word 0xaaaaaaaa, 0x0000000a, 0xcccccccc, 0xdddddddd

.align 4
result_buf:  .space 256
witness_buf: .space 256

