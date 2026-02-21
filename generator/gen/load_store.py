from __future__ import annotations

"""Generate tests for vector load/store instructions:
vle8/16/32/64, vse8/16/32/64, vlm, vsm,
vlse8/16/32/64, vsse8/16/32/64,
vle8ff/16ff/32ff/64ff,
vl1re8..vl8re64, vs1r..vs8r.
"""

from pathlib import Path

from ..common import SEWS, NUM_ELEMS, M, format_data_line, format_hex
from ..testfile import TestFile
from ..emit import VREG_DST, VREG_SRC2


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
