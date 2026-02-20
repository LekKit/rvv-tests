from __future__ import annotations

"""Generate tests for add-with-carry / subtract-with-borrow instructions:
vadc.vvm, vadc.vxm, vadc.vim,
vsbc.vvm, vsbc.vxm,
vmadc.vvm, vmadc.vxm, vmadc.vim, vmadc.vv, vmadc.vx, vmadc.vi,
vmsbc.vvm, vmsbc.vxm, vmsbc.vv, vmsbc.vx.
"""

from pathlib import Path

from ..common import SEWS, NUM_ELEMS, M, U, format_data_line, format_hex
from ..testfile import TestFile
from ..emit import VREG_DST, VREG_SRC2, VREG_SRC1
from ..compute.integer import (
    adc, sbc, madc, madc_no_carry, msbc, msbc_no_borrow,
)


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "int_adc"

    # Test data: (name, vs2[4], vs1[4], carry_mask_4bits)
    def _test_data(sew: int) -> list[tuple[str, list[int], list[int], int]]:
        m = M(sew)
        return [
            ("basic",    [1, 2, 3, 4],     [10, 20, 30, 40],  0b0000),
            ("carry",    [1, 2, 3, 4],     [10, 20, 30, 40],  0b1111),
            ("overflow", [m, m, m - 1, 0], [1, 2, 3, m],      0b1010),
            ("mixed",    [0, m, 1, m >> 1],[m, 0, m - 1, 1],  0b0101),
        ]

    # --- vadc.vvm ---
    mnemonic = "vadc.vvm"
    fname = "vadc_vvm.S"
    tf = TestFile(mnemonic, "Add with carry (vvm)")
    for sew in SEWS:
        for name, s2, s1, cmask in _test_data(sew):
            exp = []
            for i in range(NUM_ELEMS):
                cin = (cmask >> i) & 1
                exp.append(adc(s2[i], s1[i], cin, sew))
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vadc.vvm e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vadc.vvm {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {cmask}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vsbc.vvm ---
    mnemonic = "vsbc.vvm"
    fname = "vsbc_vvm.S"
    tf = TestFile(mnemonic, "Subtract with borrow (vvm)")
    for sew in SEWS:
        for name, s2, s1, bmask in _test_data(sew):
            exp = []
            for i in range(NUM_ELEMS):
                bin_ = (bmask >> i) & 1
                exp.append(sbc(s2[i], s1[i], bin_, sew))
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vsbc.vvm e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vsbc.vvm {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {bmask}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmadc.vvm (carry-out with carry-in) ---
    mnemonic = "vmadc.vvm"
    fname = "vmadc_vvm.S"
    tf = TestFile(mnemonic, "Carry-out detection with carry-in")
    for sew in SEWS:
        for name, s2, s1, cmask in _test_data(sew):
            exp_mask = 0
            for i in range(NUM_ELEMS):
                cin = (cmask >> i) & 1
                if madc(s2[i], s1[i], cin, sew):
                    exp_mask |= 1 << i
            cn = tf.next_check(f"vmadc.vvm e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vmadc.vvm {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vsm.v {VREG_DST}, (t1)")
            tf.code(f"lbu t2, 0(t1)")
            tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
            tf.code(f"li t3, {exp_mask}")
            tf.code("FAIL_IF_NE t2, t3")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {cmask}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmadc.vv (carry-out, no carry-in) ---
    mnemonic = "vmadc.vv"
    fname = "vmadc_vv.S"
    tf = TestFile(mnemonic, "Carry-out detection (no carry-in)")
    for sew in SEWS:
        m = M(sew)
        cases = [
            ("basic",    [1, 2, 3, 4],     [10, 20, 30, 40]),
            ("overflow", [m, m, m - 1, 0], [1, 2, 3, m]),
            ("max",      [m, m, m, m],     [m, m, m, m]),
        ]
        for name, s2, s1 in cases:
            exp_mask = 0
            for i in range(NUM_ELEMS):
                if madc_no_carry(s2[i], s1[i], sew):
                    exp_mask |= 1 << i
            cn = tf.next_check(f"vmadc.vv e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vmadc.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vsm.v {VREG_DST}, (t1)")
            tf.code(f"lbu t2, 0(t1)")
            tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
            tf.code(f"li t3, {exp_mask}")
            tf.code("FAIL_IF_NE t2, t3")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmsbc.vvm (borrow-out with borrow-in) ---
    mnemonic = "vmsbc.vvm"
    fname = "vmsbc_vvm.S"
    tf = TestFile(mnemonic, "Borrow-out detection with borrow-in")
    for sew in SEWS:
        for name, s2, s1, bmask in _test_data(sew):
            exp_mask = 0
            for i in range(NUM_ELEMS):
                bin_ = (bmask >> i) & 1
                if msbc(s2[i], s1[i], bin_, sew):
                    exp_mask |= 1 << i
            cn = tf.next_check(f"vmsbc.vvm e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vmsbc.vvm {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vsm.v {VREG_DST}, (t1)")
            tf.code(f"lbu t2, 0(t1)")
            tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
            tf.code(f"li t3, {exp_mask}")
            tf.code("FAIL_IF_NE t2, t3")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {bmask}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmsbc.vv (borrow-out, no borrow-in) ---
    mnemonic = "vmsbc.vv"
    fname = "vmsbc_vv.S"
    tf = TestFile(mnemonic, "Borrow-out detection (no borrow-in)")
    for sew in SEWS:
        m = M(sew)
        cases = [
            ("basic", [10, 20, 30, 40],  [1, 2, 3, 4]),
            ("under", [0, 1, 0, 2],      [1, 2, 3, 4]),
            ("equal", [5, 5, 5, 5],      [5, 5, 5, 5]),
        ]
        for name, s2, s1 in cases:
            exp_mask = 0
            for i in range(NUM_ELEMS):
                if msbc_no_borrow(s2[i], s1[i], sew):
                    exp_mask |= 1 << i
            cn = tf.next_check(f"vmsbc.vv e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_s1")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vmsbc.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vsm.v {VREG_DST}, (t1)")
            tf.code(f"lbu t2, 0(t1)")
            tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
            tf.code(f"li t3, {exp_mask}")
            tf.code("FAIL_IF_NE t2, t3")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_s1", format_data_line(s1, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vadc.vxm ---
    def _vx_test_data(sew: int) -> list[tuple[str, list[int], int, int]]:
        m = M(sew)
        return [
            ("basic",    [1, 2, 3, 4],     10,  0b0000),
            ("carry",    [1, 2, 3, 4],     10,  0b1111),
            ("overflow", [m, m, m - 1, 0], 1,   0b1010),
            ("mixed",    [0, m, 1, m >> 1],m,   0b0101),
        ]

    mnemonic = "vadc.vxm"
    fname = "vadc_vxm.S"
    tf = TestFile(mnemonic, "Add with carry (vxm)")
    for sew in SEWS:
        for name, s2, sc, cmask in _vx_test_data(sew):
            exp = []
            for i in range(NUM_ELEMS):
                cin = (cmask >> i) & 1
                exp.append(adc(s2[i], sc, cin, sew))
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vadc.vxm e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"li a0, {format_hex(sc, sew)}")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vadc.vxm {VREG_DST}, {VREG_SRC2}, a0, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {cmask}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vadc.vim ---
    mnemonic = "vadc.vim"
    fname = "vadc_vim.S"
    tf = TestFile(mnemonic, "Add with carry (vim)")
    for sew in SEWS:
        for imm in [0, 1, 5, -1, -16]:
            m = M(sew)
            s2 = [1, 2, m, 0]
            cmask = 0b1010
            exp = []
            for i in range(NUM_ELEMS):
                cin = (cmask >> i) & 1
                exp.append(adc(s2[i], imm, cin, sew))
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vadc.vim e{sew} imm={imm}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vadc.vim {VREG_DST}, {VREG_SRC2}, {imm}, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {cmask}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vsbc.vxm ---
    mnemonic = "vsbc.vxm"
    fname = "vsbc_vxm.S"
    tf = TestFile(mnemonic, "Subtract with borrow (vxm)")
    for sew in SEWS:
        for name, s2, sc, bmask in _vx_test_data(sew):
            exp = []
            for i in range(NUM_ELEMS):
                bin_ = (bmask >> i) & 1
                exp.append(sbc(s2[i], sc, bin_, sew))
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vsbc.vxm e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"li a0, {format_hex(sc, sew)}")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vsbc.vxm {VREG_DST}, {VREG_SRC2}, a0, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {bmask}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmadc.vxm (carry-out with carry-in, scalar) ---
    mnemonic = "vmadc.vxm"
    fname = "vmadc_vxm.S"
    tf = TestFile(mnemonic, "Carry-out detection with carry-in (vxm)")
    for sew in SEWS:
        for name, s2, sc, cmask in _vx_test_data(sew):
            exp_mask = 0
            for i in range(NUM_ELEMS):
                cin = (cmask >> i) & 1
                if madc(s2[i], sc, cin, sew):
                    exp_mask |= 1 << i
            cn = tf.next_check(f"vmadc.vxm e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"li a0, {format_hex(sc, sew)}")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vmadc.vxm {VREG_DST}, {VREG_SRC2}, a0, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vsm.v {VREG_DST}, (t1)")
            tf.code(f"lbu t2, 0(t1)")
            tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
            tf.code(f"li t3, {exp_mask}")
            tf.code("FAIL_IF_NE t2, t3")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {cmask}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmadc.vim (carry-out with carry-in, immediate) ---
    mnemonic = "vmadc.vim"
    fname = "vmadc_vim.S"
    tf = TestFile(mnemonic, "Carry-out detection with carry-in (vim)")
    for sew in SEWS:
        m = M(sew)
        for imm in [0, 1, 15, -1, -16]:
            s2 = [m, m - 1, 0, 1]
            cmask = 0b1010
            exp_mask = 0
            for i in range(NUM_ELEMS):
                cin = (cmask >> i) & 1
                if madc(s2[i], imm, cin, sew):
                    exp_mask |= 1 << i
            cn = tf.next_check(f"vmadc.vim e{sew} imm={imm}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vmadc.vim {VREG_DST}, {VREG_SRC2}, {imm}, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vsm.v {VREG_DST}, (t1)")
            tf.code(f"lbu t2, 0(t1)")
            tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
            tf.code(f"li t3, {exp_mask}")
            tf.code("FAIL_IF_NE t2, t3")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {cmask}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmadc.vx (carry-out, no carry-in, scalar) ---
    mnemonic = "vmadc.vx"
    fname = "vmadc_vx.S"
    tf = TestFile(mnemonic, "Carry-out detection (no carry-in, vx)")
    for sew in SEWS:
        m = M(sew)
        cases = [
            ("basic",    [1, 2, 3, 4],     10),
            ("overflow", [m, m, m - 1, 0], 1),
            ("max",      [m, m, m, m],     m),
        ]
        for name, s2, sc in cases:
            exp_mask = 0
            for i in range(NUM_ELEMS):
                if madc_no_carry(s2[i], sc, sew):
                    exp_mask |= 1 << i
            cn = tf.next_check(f"vmadc.vx e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"li a0, {format_hex(sc, sew)}")
            tf.code("SAVE_CSRS")
            tf.code(f"vmadc.vx {VREG_DST}, {VREG_SRC2}, a0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vsm.v {VREG_DST}, (t1)")
            tf.code(f"lbu t2, 0(t1)")
            tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
            tf.code(f"li t3, {exp_mask}")
            tf.code("FAIL_IF_NE t2, t3")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmadc.vi (carry-out, no carry-in, immediate) ---
    mnemonic = "vmadc.vi"
    fname = "vmadc_vi.S"
    tf = TestFile(mnemonic, "Carry-out detection (no carry-in, vi)")
    for sew in SEWS:
        m = M(sew)
        for imm in [0, 1, 15, -1]:
            s2 = [m, m - 1, 0, 1]
            exp_mask = 0
            for i in range(NUM_ELEMS):
                if madc_no_carry(s2[i], imm, sew):
                    exp_mask |= 1 << i
            cn = tf.next_check(f"vmadc.vi e{sew} imm={imm}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vmadc.vi {VREG_DST}, {VREG_SRC2}, {imm}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vsm.v {VREG_DST}, (t1)")
            tf.code(f"lbu t2, 0(t1)")
            tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
            tf.code(f"li t3, {exp_mask}")
            tf.code("FAIL_IF_NE t2, t3")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmsbc.vxm (borrow-out with borrow-in, scalar) ---
    mnemonic = "vmsbc.vxm"
    fname = "vmsbc_vxm.S"
    tf = TestFile(mnemonic, "Borrow-out detection with borrow-in (vxm)")
    for sew in SEWS:
        for name, s2, sc, bmask in _vx_test_data(sew):
            exp_mask = 0
            for i in range(NUM_ELEMS):
                bin_ = (bmask >> i) & 1
                if msbc(s2[i], sc, bin_, sew):
                    exp_mask |= 1 << i
            cn = tf.next_check(f"vmsbc.vxm e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"li a0, {format_hex(sc, sew)}")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vmsbc.vxm {VREG_DST}, {VREG_SRC2}, a0, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vsm.v {VREG_DST}, (t1)")
            tf.code(f"lbu t2, 0(t1)")
            tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
            tf.code(f"li t3, {exp_mask}")
            tf.code("FAIL_IF_NE t2, t3")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {bmask}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmsbc.vx (borrow-out, no borrow-in, scalar) ---
    mnemonic = "vmsbc.vx"
    fname = "vmsbc_vx.S"
    tf = TestFile(mnemonic, "Borrow-out detection (no borrow-in, vx)")
    for sew in SEWS:
        m = M(sew)
        cases = [
            ("basic", [10, 20, 30, 40], 1),
            ("under", [0, 1, 0, 2],     5),
            ("equal", [5, 5, 5, 5],     5),
        ]
        for name, s2, sc in cases:
            exp_mask = 0
            for i in range(NUM_ELEMS):
                if msbc_no_borrow(s2[i], sc, sew):
                    exp_mask |= 1 << i
            cn = tf.next_check(f"vmsbc.vx e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_s2")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"li a0, {format_hex(sc, sew)}")
            tf.code("SAVE_CSRS")
            tf.code(f"vmsbc.vx {VREG_DST}, {VREG_SRC2}, a0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vsm.v {VREG_DST}, (t1)")
            tf.code(f"lbu t2, 0(t1)")
            tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
            tf.code(f"li t3, {exp_mask}")
            tf.code("FAIL_IF_NE t2, t3")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(sew)
            tf.data_label(f"{tag}_s2", format_data_line(s2, sew))
    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    return generated
