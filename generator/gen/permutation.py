from __future__ import annotations

"""Generate tests for permutation instructions:
vslideup, vslidedown, vslide1up, vslide1down,
vfslide1up, vfslide1down,
vrgather, vrgatherei16, vcompress,
vmv.x.s, vmv.s.x, vfmv.f.s, vfmv.s.f,
vmv1r/2r/4r/8r, vmerge, vmv.v.
"""

from pathlib import Path

from ..common import (
    SEWS, FP_SEWS, NUM_ELEMS, M, U, format_data_line, format_hex,
    f32_to_bits, f64_to_bits,
)
from ..testfile import TestFile
from ..emit import VREG_DST, VREG_SRC2, VREG_SRC1


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "permutation"

    # --- vslideup.vi ---
    for mnemonic, direction in [("vslideup.vi", "up"), ("vslidedown.vi", "down")]:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Slide {direction} by immediate")
        for sew in SEWS:
            for offset in [0, 1, 2, 3]:
                src = [10, 20, 30, 40]
                if direction == "up":
                    # slideup: dest[i] = (i >= offset) ? src[i - offset] : dest[i]
                    # We init dest to 0, so for i < offset, dest[i] = 0
                    exp = [0] * offset + src[:NUM_ELEMS - offset]
                else:
                    # slidedown: dest[i] = (i + offset < vl) ? src[i + offset] : 0
                    exp = src[offset:] + [0] * offset
                nbytes = NUM_ELEMS * (sew // 8)
                cn = tf.next_check(f"{mnemonic} e{sew} off={offset}: result")
                tag = f"tc{cn}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code(f"la t1, {tag}_src")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                # Init dest to zero for slideup
                tf.code(f"vmv.v.i {VREG_DST}, 0")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, {offset}")
                tf.code(f"SET_TEST_NUM {cn}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(sew)
                tf.data_label(f"{tag}_src", format_data_line(src, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- vslide1up.vx / vslide1down.vx ---
    for mnemonic, direction in [("vslide1up.vx", "up"), ("vslide1down.vx", "down")]:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Slide1 {direction}")
        for sew in SEWS:
            src = [10, 20, 30, 40]
            scalar = 99
            if direction == "up":
                exp = [scalar] + src[:NUM_ELEMS - 1]
            else:
                exp = src[1:] + [scalar]
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"{mnemonic} e{sew}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"li a0, {scalar}")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, a0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(sew)
            tf.data_label(f"{tag}_src", format_data_line(src, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- vrgather.vv ---
    from ..vectors import gather_vv
    mnemonic = "vrgather.vv"
    fname = "vrgather_vv.S"
    tf = TestFile(mnemonic, "Gather elements by index vector")
    for sew in SEWS:
        for name, src, indices in gather_vv(sew):
            exp = []
            for idx in indices:
                if U(idx, sew) < NUM_ELEMS:
                    exp.append(src[U(idx, sew)])
                else:
                    exp.append(0)  # out-of-range → 0
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vrgather.vv e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vrgather.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(sew)
            tf.data_label(f"{tag}_src", format_data_line(src, sew))
            tf.data_label(f"{tag}_idx", format_data_line(indices, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

    fpath = out / fname
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vcompress.vm ---
    mnemonic = "vcompress.vm"
    fpath = out / "vcompress_vm.S"
    tf = TestFile(mnemonic, "Compress active elements")
    for sew in SEWS:
        test_cases = [
            ("all",   0b1111, [10, 20, 30, 40], [10, 20, 30, 40]),
            ("none",  0b0000, [10, 20, 30, 40], [0, 0, 0, 0]),
            ("even",  0b0101, [10, 20, 30, 40], [10, 30, 0, 0]),
            ("odd",   0b1010, [10, 20, 30, 40], [20, 40, 0, 0]),
            ("first", 0b0001, [10, 20, 30, 40], [10, 0, 0, 0]),
        ]
        for name, mask_val, src, exp in test_cases:
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vcompress e{sew} {name}: result")
            tag = f"tc{cn}"

            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            # Init dest to 0
            tf.code(f"vmv.v.i {VREG_DST}, 0")
            # Load mask into v0
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vcompress.vm {VREG_DST}, {VREG_SRC2}, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_src", format_data_line(src, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmv.x.s / vmv.s.x ---
    mnemonic = "vmv.x.s"
    fpath = out / "vmv_x_s.S"
    tf = TestFile(mnemonic, "Move vector element 0 to scalar register")
    for sew in SEWS:
        val = 0x42 & M(sew)
        cn = tf.next_check(f"vmv.x.s e{sew}: result")
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"li t2, {format_hex(val, sew)}")
        tf.code(f"vmv.v.x {VREG_SRC2}, t2")
        tf.code("SAVE_CSRS")
        tf.code(f"vmv.x.s t3, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"li t4, {format_hex(val, sew)}")
        tf.code("FAIL_IF_NE t3, t4")
        tf.code("CHECK_CSRS_UNCHANGED")
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmv.v.v / vmv.v.x / vmv.v.i ---
    mnemonic = "vmv.v.v"
    fpath = out / "vmv_v_v.S"
    tf = TestFile(mnemonic, "Vector-vector move")
    for sew in SEWS:
        src = [10, 20, 30, 40]
        nbytes = NUM_ELEMS * (sew // 8)
        cn = tf.next_check(f"vmv.v.v e{sew}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_src")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vmv.v.v {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_src, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data_align(sew)
        tf.data_label(f"{tag}_src", format_data_line(src, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmerge.vvm / vmerge.vxm / vmerge.vim ---
    mnemonic = "vmerge.vvm"
    fpath = out / "vmerge_vvm.S"
    tf = TestFile(mnemonic, "Vector merge with mask")
    for sew in SEWS:
        vs2 = [10, 20, 30, 40]
        vs1 = [100, 200, 300, 400]
        mask_val = 0b1010  # select vs1 for elements 1,3; vs2 for 0,2
        exp = []
        for i in range(NUM_ELEMS):
            if mask_val & (1 << i):
                exp.append(vs1[i])
            else:
                exp.append(vs2[i])
        nbytes = NUM_ELEMS * (sew // 8)
        cn = tf.next_check(f"vmerge.vvm e{sew}: result")
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
        tf.code(f"vmerge.vvm {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}, v0")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data(".align 1")
        tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
        tf.data_align(sew)
        tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
        tf.data_label(f"{tag}_s1", format_data_line(vs1, sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmv.s.x  (scalar insert into element 0) ---
    mnemonic = "vmv.s.x"
    fpath = out / "vmv_s_x.S"
    tf = TestFile(mnemonic, "Move scalar to vector element 0")
    for sew in SEWS:
        val = 0x42 & M(sew)
        cn = tf.next_check(f"vmv.s.x e{sew}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"vmv.v.i {VREG_DST}, 0")  # clear dest
        tf.code(f"li t2, {format_hex(val, sew)}")
        tf.code("SAVE_CSRS")
        tf.code(f"vmv.s.x {VREG_DST}, t2")
        tf.code(f"SET_TEST_NUM {cn}")
        # Read back element 0
        tf.code(f"vmv.x.s t3, {VREG_DST}")
        tf.code(f"li t4, {format_hex(val, sew)}")
        tf.code("FAIL_IF_NE t3, t4")
        tf.code("CHECK_CSRS_UNCHANGED")
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vfmv.f.s  (FP element 0 → FP register) ---
    mnemonic = "vfmv.f.s"
    fpath = out / "vfmv_f_s.S"
    tf = TestFile(mnemonic, "Move vector FP element 0 to FP register")
    for sew in FP_SEWS:
        b = f32_to_bits if sew == 32 else f64_to_bits
        val_bits = b(3.14) if sew == 32 else b(3.14)
        load_insn = "flw" if sew == 32 else "fld"
        store_insn = "fsw" if sew == 32 else "fsd"
        cn = tf.next_check(f"vfmv.f.s e{sew}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_src")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vfmv.f.s fa1, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        # Store fa1 to memory and compare
        tf.code(f"la t1, result_buf")
        tf.code(f"{store_insn} fa1, 0(t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {sew // 8}")
        tf.code("CHECK_CSRS_UNCHANGED_FP")
        src_vals = [val_bits, 0, 0, 0]
        tf.data_align(sew)
        tf.data_label(f"{tag}_src", format_data_line(src_vals, sew))
        tf.data_label(f"{tag}_exp", format_data_line([val_bits], sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vfmv.s.f  (FP register → element 0) ---
    mnemonic = "vfmv.s.f"
    fpath = out / "vfmv_s_f.S"
    tf = TestFile(mnemonic, "Move FP scalar to vector element 0")
    for sew in FP_SEWS:
        b = f32_to_bits if sew == 32 else f64_to_bits
        val_bits = b(2.718) if sew == 32 else b(2.718)
        load_insn = "flw" if sew == 32 else "fld"
        cn = tf.next_check(f"vfmv.s.f e{sew}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"vmv.v.i {VREG_DST}, 0")
        tf.code(f"la t1, {tag}_fsc")
        tf.code(f"{load_insn} fa0, 0(t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vfmv.s.f {VREG_DST}, fa0")
        tf.code(f"SET_TEST_NUM {cn}")
        # Store result and check element 0
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {sew // 8}")
        tf.code("CHECK_CSRS_UNCHANGED_FP")
        tf.data_align(sew)
        tf.data_label(f"{tag}_fsc", format_data_line([val_bits], sew))
        tf.data_label(f"{tag}_exp", format_data_line([val_bits], sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vfmerge.vfm  (FP merge with mask) ---
    mnemonic = "vfmerge.vfm"
    fpath = out / "vfmerge_vfm.S"
    tf = TestFile(mnemonic, "FP vector merge with mask")
    for sew in FP_SEWS:
        b = f32_to_bits if sew == 32 else f64_to_bits
        load_insn = "flw" if sew == 32 else "fld"
        vs2 = [b(1.0), b(2.0), b(3.0), b(4.0)]
        scalar_bits = b(99.0)
        mask_val = 0b1010
        exp = []
        for i in range(NUM_ELEMS):
            if mask_val & (1 << i):
                exp.append(scalar_bits)
            else:
                exp.append(vs2[i])
        nbytes = NUM_ELEMS * (sew // 8)
        cn = tf.next_check(f"vfmerge.vfm e{sew}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_s2")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code(f"la t1, {tag}_fsc")
        tf.code(f"{load_insn} fa0, 0(t1)")
        tf.code(f"la t1, {tag}_mask")
        tf.code(f"vlm.v v0, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vfmerge.vfm {VREG_DST}, {VREG_SRC2}, fa0, v0")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED_FP")
        tf.data(".align 1")
        tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
        tf.data_align(sew)
        tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
        tf.data_label(f"{tag}_fsc", format_data_line([scalar_bits], sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vfmv.v.f  (FP scalar broadcast) ---
    mnemonic = "vfmv.v.f"
    fpath = out / "vfmv_v_f.S"
    tf = TestFile(mnemonic, "FP scalar broadcast to vector")
    for sew in FP_SEWS:
        b = f32_to_bits if sew == 32 else f64_to_bits
        load_insn = "flw" if sew == 32 else "fld"
        val_bits = b(7.5)
        exp = [val_bits] * NUM_ELEMS
        nbytes = NUM_ELEMS * (sew // 8)
        cn = tf.next_check(f"vfmv.v.f e{sew}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_fsc")
        tf.code(f"{load_insn} fa0, 0(t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vfmv.v.f {VREG_DST}, fa0")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED_FP")
        tf.data_align(sew)
        tf.data_label(f"{tag}_fsc", format_data_line([val_bits], sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vrgather.vx ---
    mnemonic = "vrgather.vx"
    fpath = out / "vrgather_vx.S"
    tf = TestFile(mnemonic, "Gather by scalar index")
    for sew in SEWS:
        src = [10, 20, 30, 40]
        for idx in range(NUM_ELEMS):
            exp = [src[idx]] * NUM_ELEMS  # all elements get src[idx]
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vrgather.vx e{sew} idx={idx}: result")
            tag = f"tc{cn}"
            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"li a0, {idx}")
            tf.code("SAVE_CSRS")
            tf.code(f"vrgather.vx {VREG_DST}, {VREG_SRC2}, a0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data_align(sew)
            tf.data_label(f"{tag}_src", format_data_line(src, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vrgather.vi ---
    mnemonic = "vrgather.vi"
    fpath = out / "vrgather_vi.S"
    tf = TestFile(mnemonic, "Gather by immediate index")
    for sew in SEWS:
        src = [10, 20, 30, 40]
        for idx in range(NUM_ELEMS):
            exp = [src[idx]] * NUM_ELEMS
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vrgather.vi e{sew} idx={idx}: result")
            tag = f"tc{cn}"
            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"vrgather.vi {VREG_DST}, {VREG_SRC2}, {idx}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data_align(sew)
            tf.data_label(f"{tag}_src", format_data_line(src, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vrgatherei16.vv ---
    mnemonic = "vrgatherei16.vv"
    fpath = out / "vrgatherei16_vv.S"
    tf = TestFile(mnemonic, "Gather with 16-bit index vector")
    for sew in SEWS:
        src = [10, 20, 30, 40]
        test_cases_g16 = [
            ("identity", [0, 1, 2, 3]),
            ("reverse",  [3, 2, 1, 0]),
            ("splat",    [0, 0, 0, 0]),
        ]
        for name, indices in test_cases_g16:
            exp = []
            for idx in indices:
                exp.append(src[idx] if idx < NUM_ELEMS else 0)
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vrgatherei16 e{sew} {name}: result")
            tag = f"tc{cn}"
            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            # Load 16-bit indices
            tf.code(f"vsetvli t0, t0, e16, m1, tu, mu")
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle16.v {VREG_SRC1}, (t1)")
            # Restore vtype for the operation
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code("SAVE_CSRS")
            tf.code(f"vrgatherei16.vv {VREG_DST}, {VREG_SRC2}, {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data_align(sew)
            tf.data_label(f"{tag}_src", format_data_line(src, sew))
            tf.data(f".align 1")
            tf.data_label(f"{tag}_idx", format_data_line(indices, 16))
            tf.data_align(sew)
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vfslide1up.vf / vfslide1down.vf ---
    for mnemonic, direction in [("vfslide1up.vf", "up"),
                                ("vfslide1down.vf", "down")]:
        fname = mnemonic.replace(".", "_") + ".S"
        fpath = out / fname
        tf = TestFile(mnemonic, f"FP slide1 {direction}")
        for sew in FP_SEWS:
            b = f32_to_bits if sew == 32 else f64_to_bits
            load_insn = "flw" if sew == 32 else "fld"
            src = [b(1.0), b(2.0), b(3.0), b(4.0)]
            scalar_bits = b(99.0)
            if direction == "up":
                exp = [scalar_bits] + src[:NUM_ELEMS - 1]
            else:
                exp = src[1:] + [scalar_bits]
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"{mnemonic} e{sew}: result")
            tag = f"tc{cn}"
            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code(f"la t1, {tag}_fsc")
            tf.code(f"{load_insn} fa0, 0(t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, fa0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED_FP")
            tf.data_align(sew)
            tf.data_label(f"{tag}_src", format_data_line(src, sew))
            tf.data_label(f"{tag}_fsc", format_data_line([scalar_bits], sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
        tf.write(fpath)
        generated.append(str(fpath))

    # --- vslideup.vx / vslidedown.vx ---
    for mnemonic, direction in [("vslideup.vx", "up"), ("vslidedown.vx", "down")]:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Slide {direction} by scalar register")
        for sew in SEWS:
            src = [10, 20, 30, 40]
            for offset in [0, 1, 2, 3]:
                if direction == "up":
                    exp = [0] * offset + src[:NUM_ELEMS - offset]
                else:
                    exp = src[offset:] + [0] * offset
                nbytes = NUM_ELEMS * (sew // 8)
                cn = tf.next_check(f"{mnemonic} e{sew} off={offset}: result")
                tag = f"tc{cn}"

                tf.blank()
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
                tf.code(f"la t1, {tag}_src")
                tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
                tf.code(f"vmv.v.i {VREG_DST}, 0")
                tf.code(f"li a0, {offset}")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}, a0")
                tf.code(f"SET_TEST_NUM {cn}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code("CHECK_CSRS_UNCHANGED")

                tf.data_align(sew)
                tf.data_label(f"{tag}_src", format_data_line(src, sew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, sew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- vmerge.vxm ---
    mnemonic = "vmerge.vxm"
    fpath = out / "vmerge_vxm.S"
    tf = TestFile(mnemonic, "Vector merge with scalar and mask")
    for sew in SEWS:
        vs2 = [10, 20, 30, 40]
        scalar = 99
        mask_val = 0b1010
        exp = []
        for i in range(NUM_ELEMS):
            if mask_val & (1 << i):
                exp.append(scalar)
            else:
                exp.append(vs2[i])
        nbytes = NUM_ELEMS * (sew // 8)
        cn = tf.next_check(f"vmerge.vxm e{sew}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_s2")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code(f"li a0, {scalar}")
        tf.code(f"la t1, {tag}_mask")
        tf.code(f"vlm.v v0, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vmerge.vxm {VREG_DST}, {VREG_SRC2}, a0, v0")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data(".align 1")
        tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
        tf.data_align(sew)
        tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmerge.vim ---
    mnemonic = "vmerge.vim"
    fpath = out / "vmerge_vim.S"
    tf = TestFile(mnemonic, "Vector merge with immediate and mask")
    for sew in SEWS:
        vs2 = [10, 20, 30, 40]
        imm = 5
        mask_val = 0b0101
        exp = []
        for i in range(NUM_ELEMS):
            if mask_val & (1 << i):
                exp.append(imm)
            else:
                exp.append(vs2[i])
        nbytes = NUM_ELEMS * (sew // 8)
        cn = tf.next_check(f"vmerge.vim e{sew}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_s2")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code(f"la t1, {tag}_mask")
        tf.code(f"vlm.v v0, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vmerge.vim {VREG_DST}, {VREG_SRC2}, {imm}, v0")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data(".align 1")
        tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
        tf.data_align(sew)
        tf.data_label(f"{tag}_s2", format_data_line(vs2, sew))
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmv.v.x  (scalar broadcast) ---
    mnemonic = "vmv.v.x"
    fpath = out / "vmv_v_x.S"
    tf = TestFile(mnemonic, "Scalar broadcast to vector")
    for sew in SEWS:
        val = 42 & M(sew)
        exp = [val] * NUM_ELEMS
        nbytes = NUM_ELEMS * (sew // 8)
        cn = tf.next_check(f"vmv.v.x e{sew}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"li a0, {format_hex(val, sew)}")
        tf.code("SAVE_CSRS")
        tf.code(f"vmv.v.x {VREG_DST}, a0")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data_align(sew)
        tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmv.v.i  (immediate broadcast) ---
    mnemonic = "vmv.v.i"
    fpath = out / "vmv_v_i.S"
    tf = TestFile(mnemonic, "Immediate broadcast to vector")
    for sew in SEWS:
        for imm in [0, 1, 5, -1, 15, -16]:
            # simm5 is sign-extended to sew bits
            val = imm & M(sew)
            exp = [val] * NUM_ELEMS
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"vmv.v.i e{sew} imm={imm}: result")
            tag = f"tc{cn}"
            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code("SAVE_CSRS")
            tf.code(f"vmv.v.i {VREG_DST}, {imm}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data_align(sew)
            tf.data_label(f"{tag}_exp", format_data_line(exp, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- Whole register move: vmv1r.v ---
    mnemonic = "vmv1r.v"
    fpath = out / "vmv1r_v.S"
    tf = TestFile(mnemonic, "Whole register move (1 register)")
    sew = 32  # Whole-reg moves are vtype-independent, use e32 for verification
    src = [0xDEADBEEF, 0xCAFEBABE, 0x12345678, 0x9ABCDEF0]
    nbytes = NUM_ELEMS * (sew // 8)
    cn = tf.next_check("vmv1r.v: result")
    tag = f"tc{cn}"
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
    tf.code(f"la t1, {tag}_src")
    tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
    tf.code("SAVE_CSRS")
    tf.code(f"vmv1r.v {VREG_DST}, {VREG_SRC2}")
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
    tf.code(f"CHECK_MEM result_buf, {tag}_src, {nbytes}")
    tf.code("CHECK_CSRS_UNCHANGED")
    tf.data_align(sew)
    tf.data_label(f"{tag}_src", format_data_line(src, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- Whole register move: vmv2r.v / vmv4r.v / vmv8r.v ---
    # These copy NR consecutive vector registers.
    # Register allocation:
    #   vmv2r: src=v16(v16-v17), dst=v8(v8-v9)   → LMUL=m2
    #   vmv4r: src=v16(v16-v19), dst=v8(v8-v11)   → LMUL=m4
    #   vmv8r: src=v16(v16-v23), dst=v8(v8-v15)   → LMUL=m8
    # We use e32 for loading/storing, so each register holds NUM_ELEMS
    # (4) e32 elements.  Total elements checked = NUM_ELEMS * NR.
    for nr, lmul_str in [(2, "m2"), (4, "m4"), (8, "m8")]:
        mnemonic = f"vmv{nr}r.v"
        fpath = out / f"vmv{nr}r_v.S"
        tf = TestFile(mnemonic, f"Whole register move ({nr} registers)")
        sew = 32
        # Build test data: NR * NUM_ELEMS distinct e32 values
        total_elems = NUM_ELEMS * nr
        src_data = [(0xA0000000 + i * 0x11111111) & M(sew) for i in range(total_elems)]
        nbytes = total_elems * (sew // 8)
        cn = tf.next_check(f"vmv{nr}r.v: result")
        tag = f"tc{cn}"

        # Set LMUL to NR so we can load all NR regs at once
        tf.code(f"li t0, {total_elems}")
        tf.code(f"vsetvli t0, t0, e{sew}, {lmul_str}, tu, mu")
        # Load source
        tf.code(f"la t1, {tag}_src")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        # Store result
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_src, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")

        tf.data_align(sew)
        tf.data_label(f"{tag}_src", format_data_line(src_data, sew))
        tf.write(fpath)
        generated.append(str(fpath))

    return generated
