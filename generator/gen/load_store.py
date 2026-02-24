from __future__ import annotations

"""Generate tests for vector load/store instructions:
vle8/16/32/64, vse8/16/32/64, vlm, vsm,
vlse8/16/32/64, vsse8/16/32/64,
vle8ff/16ff/32ff/64ff,
vl1re8..vl8re64, vs1r..vs8r.
"""

from pathlib import Path

from ..common import SEWS, NUM_ELEMS, M, format_data_line, format_hex, witness_pattern
from ..testfile import TestFile
from ..emit import VREG_DST, VREG_SRC2, VREG_WITNESS


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []

    # --- Unit-stride loads + stores ---
    for sew in SEWS:
        out_load = base_dir / "tests" / "load"
        out_store = base_dir / "tests" / "store"

        # vle{sew}.v
        mnemonic = f"vle{sew}.v"
        fpath = out_load / f"vle{sew}.S"
        tf = TestFile(mnemonic, f"Unit-stride load of {sew}-bit elements")

        src_vals = [i + 1 for i in range(NUM_ELEMS)]
        nbytes = NUM_ELEMS * (sew // 8)

        # Test: load data, store to different buffer, compare with original
        cn = tf.next_check(f"vle{sew}.v: basic load-store round-trip")
        tag = f"tc{cn}"
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_data")
        tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_data, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")

        # Test: verify element values using scalar loads
        cn2 = tf.next_check(f"vle{sew}.v: scalar verify element 0")
        tf.blank()
        tf.code(f"SET_TEST_NUM {cn2}")
        scalar_load = {8: "lbu", 16: "lhu", 32: "lwu", 64: "ld"}[sew]
        tf.code(f"la t1, result_buf")
        tf.code(f"{scalar_load} t2, 0(t1)")
        tf.code(f"li t3, {src_vals[0]}")
        tf.code("FAIL_IF_NE t2, t3")

        tf.data_align(sew)
        tf.data_label(f"{tag}_data", format_data_line(src_vals, sew))

        # Masked load: mask=0b0101, elements 0,2 loaded, 1,3 preserved from vd_init
        mask_bits = 0b0101
        vd_init = [0xDD & M(sew)] * NUM_ELEMS
        mem_data = [0xAA, 0xBB, 0xCC, 0xDD]
        mem_data = [v & M(sew) for v in mem_data]
        exp_masked = []
        for i in range(NUM_ELEMS):
            if mask_bits & (1 << i):
                exp_masked.append(mem_data[i])
            else:
                exp_masked.append(vd_init[i])
        cn_m = tf.next_check(f"vle{sew}.v masked: active+inactive elements")
        tag_m = f"tc{cn_m}"
        tf.blank()
        tf.comment(f"Masked load: mask=0b{mask_bits:04b}, tu/mu policy")
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag_m}_vd")
        tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
        tf.code(f"la t1, {tag_m}_mask")
        tf.code("vlm.v v0, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"la t1, {tag_m}_mem")
        tf.code(f"vle{sew}.v {VREG_DST}, (t1), v0.t")
        tf.code(f"SET_TEST_NUM {cn_m}")
        tf.code("la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag_m}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")

        tf.data(".align 1")
        tf.data_label(f"{tag_m}_mask", f"    .byte {mask_bits}")
        tf.data_align(sew)
        tf.data_label(f"{tag_m}_vd", format_data_line(vd_init, sew))
        tf.data_label(f"{tag_m}_mem", format_data_line(mem_data, sew))
        tf.data_label(f"{tag_m}_exp", format_data_line(exp_masked, sew))

        tf.write(fpath)
        generated.append(str(fpath))

        # vse{sew}.v  (store test)
        mnemonic = f"vse{sew}.v"
        fpath = out_store / f"vse{sew}.S"
        tf = TestFile(mnemonic, f"Unit-stride store of {sew}-bit elements")

        cn = tf.next_check(f"vse{sew}.v: store and scalar-verify")
        tag = f"tc{cn}"
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        # Load known data
        tf.code(f"la t1, {tag}_data")
        tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
        # Store to buffer
        tf.code("SAVE_CSRS")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"SET_TEST_NUM {cn}")
        # Verify each byte with scalar
        scalar_load = {8: "lbu", 16: "lhu", 32: "lwu", 64: "ld"}[sew]
        elem_bytes = sew // 8
        for i, val in enumerate(src_vals):
            tf.code(f"{scalar_load} t2, {i * elem_bytes}(t1)")
            tf.code(f"li t3, {val}")
            tf.code("FAIL_IF_NE t2, t3")
        tf.code("CHECK_CSRS_UNCHANGED")

        tf.data_align(sew)
        tf.data_label(f"{tag}_data", format_data_line(src_vals, sew))

        # Masked store: mask=0b0101, only elements 0,2 written to buffer
        mask_bits = 0b0101
        store_data = [0x11, 0x22, 0x33, 0x44]
        store_data = [v & M(sew) for v in store_data]
        buf_init = [0xFF & M(sew)] * NUM_ELEMS
        exp_store = []
        for i in range(NUM_ELEMS):
            if mask_bits & (1 << i):
                exp_store.append(store_data[i])
            else:
                exp_store.append(buf_init[i])
        cn_m = tf.next_check(f"vse{sew}.v masked: only active elements stored")
        tag_m = f"tc{cn_m}"
        tf.blank()
        tf.comment(f"Masked store: mask=0b{mask_bits:04b}, only active written")
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag_m}_src")
        tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
        tf.code(f"la t1, {tag_m}_mask")
        tf.code("vlm.v v0, (t1)")
        # Pre-fill result_buf with buf_init pattern
        tf.code(f"la t1, {tag_m}_buf_init")
        tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
        tf.code("la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_SRC2}, (t1)")
        # Now do masked store
        tf.code("SAVE_CSRS")
        tf.code("la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1), v0.t")
        tf.code(f"SET_TEST_NUM {cn_m}")
        tf.code(f"CHECK_MEM result_buf, {tag_m}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")

        tf.data(".align 1")
        tf.data_label(f"{tag_m}_mask", f"    .byte {mask_bits}")
        tf.data_align(sew)
        tf.data_label(f"{tag_m}_src", format_data_line(store_data, sew))
        tf.data_label(f"{tag_m}_buf_init", format_data_line(buf_init, sew))
        tf.data_label(f"{tag_m}_exp", format_data_line(exp_store, sew))

        tf.write(fpath)
        generated.append(str(fpath))

    # --- vlm.v / vsm.v (mask load/store) ---
    fpath = base_dir / "tests" / "load" / "vlm.S"
    tf = TestFile("vlm.v", "Unit-stride mask load")
    cn = tf.next_check("vlm.v: basic mask load")
    tag = f"tc{cn}"
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e8, m1, tu, mu")
    tf.code(f"la t1, {tag}_mask")
    tf.code(f"vlm.v {VREG_DST}, (t1)")
    tf.code("SAVE_CSRS")
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code(f"la t1, result_buf")
    tf.code(f"vsm.v {VREG_DST}, (t1)")
    tf.code(f"lbu t2, 0(t1)")
    tf.code(f"andi t2, t2, 0xF")
    tf.code(f"li t3, 0b1010")
    tf.code("FAIL_IF_NE t2, t3")
    tf.code("CHECK_CSRS_UNCHANGED")
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", "    .byte 0b1010")
    tf.write(fpath)
    generated.append(str(fpath))

    # --- Strided stores ---
    for sew in SEWS:
        mnemonic = f"vsse{sew}.v"
        fpath = base_dir / "tests" / "store" / f"vsse{sew}.S"
        tf = TestFile(mnemonic, f"Strided store of {sew}-bit elements")
        elem_bytes = sew // 8
        stride = elem_bytes * 2
        src_vals = [10, 20, 30, 40]
        # After strided store with stride=2*elem_bytes, data at offsets
        # 0, 2*eb, 4*eb, 6*eb should contain src_vals
        cn = tf.next_check(f"vsse{sew}.v stride={stride}: result")
        tag = f"tc{cn}"
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_data")
        tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
        # Clear result buffer (8 elements worth)
        tf.code(f"la t1, result_buf")
        tf.code(f"li t2, {stride}")
        tf.code("SAVE_CSRS")
        tf.code(f"vsse{sew}.v {VREG_DST}, (t1), t2")
        tf.code(f"SET_TEST_NUM {cn}")
        # Verify by loading with same stride
        tf.code(f"vlse{sew}.v {VREG_SRC2}, (t1), t2")
        tf.code(f"la t1, witness_buf")
        tf.code(f"vse{sew}.v {VREG_SRC2}, (t1)")
        nbytes = NUM_ELEMS * elem_bytes
        tf.code(f"CHECK_MEM witness_buf, {tag}_data, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data_align(sew)
        tf.data_label(f"{tag}_data", format_data_line(src_vals, sew))

        # Masked strided store: mask=0b0101, only elements 0,2 stored
        mask_bits = 0b0101
        store_data_m = [0x11, 0x22, 0x33, 0x44]
        store_data_m = [v & M(sew) for v in store_data_m]
        cn_m = tf.next_check(f"vsse{sew}.v stride={stride} masked: result")
        tag_m = f"tc{cn_m}"
        tf.blank()
        tf.comment(f"Masked strided store: mask=0b{mask_bits:04b}")
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag_m}_src")
        tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
        tf.code(f"la t1, {tag_m}_mask")
        tf.code("vlm.v v0, (t1)")
        # Clear result_buf (zero fill via vmv.v.i + strided store)
        tf.code(f"vmv.v.i {VREG_SRC2}, 0")
        tf.code("la t1, result_buf")
        tf.code(f"li t2, {stride}")
        tf.code(f"vsse{sew}.v {VREG_SRC2}, (t1), t2")
        # Now masked strided store
        tf.code("SAVE_CSRS")
        tf.code("la t1, result_buf")
        tf.code(f"vsse{sew}.v {VREG_DST}, (t1), t2, v0.t")
        tf.code(f"SET_TEST_NUM {cn_m}")
        # Verify: load back with stride and check
        tf.code(f"vlse{sew}.v {VREG_SRC2}, (t1), t2")
        tf.code("la t1, witness_buf")
        tf.code(f"vse{sew}.v {VREG_SRC2}, (t1)")
        # Expected: elements 0,2 = store_data_m, elements 1,3 = 0
        exp_strided_m = []
        for i in range(NUM_ELEMS):
            if mask_bits & (1 << i):
                exp_strided_m.append(store_data_m[i])
            else:
                exp_strided_m.append(0)
        tf.code(f"CHECK_MEM witness_buf, {tag_m}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data(".align 1")
        tf.data_label(f"{tag_m}_mask", f"    .byte {mask_bits}")
        tf.data_align(sew)
        tf.data_label(f"{tag_m}_src", format_data_line(store_data_m, sew))
        tf.data_label(f"{tag_m}_exp", format_data_line(exp_strided_m, sew))

        tf.write(fpath)
        generated.append(str(fpath))

    # --- Strided loads ---
    for sew in SEWS:
        mnemonic = f"vlse{sew}.v"
        fpath = base_dir / "tests" / "load" / f"vlse{sew}.S"
        tf = TestFile(mnemonic, f"Strided load of {sew}-bit elements")

        elem_bytes = sew // 8
        stride = elem_bytes * 2  # stride = 2 elements apart
        # Source data: elements at offsets 0, 2*eb, 4*eb, 6*eb
        # We need at least stride * (NUM_ELEMS-1) + elem_bytes bytes of source
        # Use 8 elements worth of data, pick every other one
        full_data = list(range(1, 9))
        expected = [full_data[i * 2] for i in range(NUM_ELEMS)]
        nbytes = NUM_ELEMS * elem_bytes

        cn = tf.next_check(f"vlse{sew}.v stride={stride}: result")
        tag = f"tc{cn}"
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_data")
        tf.code(f"li t2, {stride}")
        tf.code("SAVE_CSRS")
        tf.code(f"vlse{sew}.v {VREG_DST}, (t1), t2")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")

        tf.data_align(sew)
        tf.data_label(f"{tag}_data", format_data_line(full_data, sew))
        tf.data_label(f"{tag}_exp", format_data_line(expected, sew))

        # Masked strided load: mask=0b0101, elements 0,2 loaded, 1,3 preserved
        mask_bits = 0b0101
        vd_init_m = [0xDD & M(sew)] * NUM_ELEMS
        exp_masked_sl = []
        for i in range(NUM_ELEMS):
            if mask_bits & (1 << i):
                exp_masked_sl.append(full_data[i * 2])
            else:
                exp_masked_sl.append(vd_init_m[i])
        cn_m = tf.next_check(f"vlse{sew}.v stride={stride} masked: result")
        tag_m = f"tc{cn_m}"
        tf.blank()
        tf.comment(f"Masked strided load: mask=0b{mask_bits:04b}")
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag_m}_vd")
        tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
        tf.code(f"la t1, {tag_m}_mask")
        tf.code("vlm.v v0, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"la t1, {tag_m}_data")
        tf.code(f"li t2, {stride}")
        tf.code(f"vlse{sew}.v {VREG_DST}, (t1), t2, v0.t")
        tf.code(f"SET_TEST_NUM {cn_m}")
        tf.code("la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag_m}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data(".align 1")
        tf.data_label(f"{tag_m}_mask", f"    .byte {mask_bits}")
        tf.data_align(sew)
        tf.data_label(f"{tag_m}_vd", format_data_line(vd_init_m, sew))
        tf.data_label(f"{tag_m}_data", format_data_line(full_data, sew))
        tf.data_label(f"{tag_m}_exp", format_data_line(exp_masked_sl, sew))

        tf.write(fpath)
        generated.append(str(fpath))

    # --- Fault-only-first loads ---
    for sew in SEWS:
        mnemonic = f"vle{sew}ff.v"
        fpath = base_dir / "tests" / "load" / f"vle{sew}ff.S"
        tf = TestFile(mnemonic, f"Fault-only-first load {sew}-bit")

        src_vals = [10, 20, 30, 40]
        nbytes = NUM_ELEMS * (sew // 8)

        # Basic test: should behave like normal load when no faults
        cn = tf.next_check(f"vle{sew}ff.v: basic (no fault)")
        tag = f"tc{cn}"
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_data")
        tf.code("SAVE_CSRS")
        tf.code(f"vle{sew}ff.v {VREG_DST}, (t1)")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_data, {nbytes}")
        # vl may have been reduced by ff, but for valid memory it shouldn't
        tf.code("CHECK_CSRS_UNCHANGED")

        tf.data_align(sew)
        tf.data_label(f"{tag}_data", format_data_line(src_vals, sew))

        tf.write(fpath)
        generated.append(str(fpath))

    # --- Whole register load/store ---
    for nregs in [1, 2, 4, 8]:
        for sew in SEWS:
            mnemonic = f"vl{nregs}re{sew}.v"
            fpath = base_dir / "tests" / "load" / f"vl{nregs}re{sew}.S"
            tf = TestFile(mnemonic, f"Whole register load ({nregs} reg, eew={sew})")

            # Use e32/m1 for setup, load nregs worth of data
            src_vals = list(range(1, NUM_ELEMS + 1))
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"{mnemonic}: round-trip")
            tag = f"tc{cn}"

            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_data")
            tf.code("SAVE_CSRS")
            tf.code(f"vl{nregs}re{sew}.v {VREG_DST}, (t1)")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            # Only check first register's worth of data
            tf.code(f"CHECK_MEM result_buf, {tag}_data, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(sew)
            tf.data_label(f"{tag}_data", format_data_line(src_vals, sew))

            tf.write(fpath)
            generated.append(str(fpath))

    # --- Whole register store ---
    for nregs in [1, 2, 4, 8]:
        mnemonic = f"vs{nregs}r.v"
        fpath = base_dir / "tests" / "store" / f"vs{nregs}r.S"
        tf = TestFile(mnemonic, f"Whole register store ({nregs} reg)")
        sew = 32
        src_vals = [0xDEADBEEF, 0xCAFEBABE, 0x12345678, 0x9ABCDEF0]
        nbytes = NUM_ELEMS * (sew // 8)
        cn = tf.next_check(f"vs{nregs}r.v: round-trip")
        tag = f"tc{cn}"

        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code(f"la t1, {tag}_data")
        tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"la t1, result_buf")
        tf.code(f"vs{nregs}r.v {VREG_DST}, (t1)")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"CHECK_MEM result_buf, {tag}_data, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")

        tf.data_align(sew)
        tf.data_label(f"{tag}_data", format_data_line(src_vals, sew))

        tf.write(fpath)
        generated.append(str(fpath))

    # --- vsm.v (standalone mask store test) ---
    fpath = base_dir / "tests" / "store" / "vsm.S"
    tf = TestFile("vsm.v", "Unit-stride mask store")
    cn = tf.next_check("vsm.v: basic mask store")
    tag = f"tc{cn}"
    tf.code(f"li t0, {NUM_ELEMS}")
    tf.code(f"vsetvli t0, t0, e8, m1, tu, mu")
    # Load a known mask pattern into a vector register
    tf.code(f"la t1, {tag}_mask")
    tf.code(f"vlm.v {VREG_DST}, (t1)")
    tf.code("SAVE_CSRS")
    # Store it using vsm.v
    tf.code(f"la t1, result_buf")
    tf.code(f"vsm.v {VREG_DST}, (t1)")
    tf.code(f"SET_TEST_NUM {cn}")
    # Verify by scalar load
    tf.code(f"lbu t2, 0(t1)")
    tf.code(f"andi t2, t2, 0xF")  # only low NUM_ELEMS bits
    tf.code(f"li t3, 0b0101")
    tf.code("FAIL_IF_NE t2, t3")
    tf.code("CHECK_CSRS_UNCHANGED")
    tf.data(".align 1")
    tf.data_label(f"{tag}_mask", "    .byte 0b0101")
    tf.write(fpath)
    generated.append(str(fpath))

    return generated
