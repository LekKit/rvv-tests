from __future__ import annotations

"""Generate tests for indexed (vector-offset) load/store instructions:
vluxei8/16/32/64, vloxei8/16/32/64 (indexed loads),
vsuxei8/16/32/64, vsoxei8/16/32/64 (indexed stores).

vluxei = unordered indexed load, vloxei = ordered indexed load.
Functionally identical for correctness testing; only ordering differs.
vsuxei = unordered indexed store, vsoxei = ordered indexed store.
"""

from pathlib import Path

from ..common import SEWS, NUM_ELEMS, M, format_data_line, format_hex
from ..testfile import TestFile
from ..emit import VREG_DST, VREG_SRC2, VREG_SRC1


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []

    # --- Indexed loads: vluxei / vloxei ---
    # Index EEW (the width of the index elements) can differ from data SEW.
    # For simplicity, we test index_eew == data_sew (matching widths).
    # The index vector contains byte offsets into the source buffer.
    for prefix, subdir_name in [("vluxei", "load"), ("vloxei", "load")]:
        for sew in SEWS:
            mnemonic = f"{prefix}{sew}.v"
            fname = f"{prefix}{sew}.S"
            fpath = base_dir / "tests" / subdir_name / fname
            tf = TestFile(mnemonic, f"Indexed load {mnemonic}")

            elem_bytes = sew // 8
            # Source data: 8 elements (extra for out-of-order access)
            src_data = list(range(10, 10 + 8))

            # Test 1: identity (sequential) access
            indices_identity = [i * elem_bytes for i in range(NUM_ELEMS)]
            exp_identity = src_data[:NUM_ELEMS]
            cn = tf.next_check(f"{mnemonic} identity: result")
            tag = f"tc{cn}"
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_data")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, (t1), {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            nbytes = NUM_ELEMS * elem_bytes
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data_align(sew)
            tf.data_label(f"{tag}_data", format_data_line(src_data, sew))
            tf.data_label(f"{tag}_idx", format_data_line(indices_identity, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp_identity, sew))

            # Test 2: reverse access
            indices_reverse = [
                (NUM_ELEMS - 1 - i) * elem_bytes for i in range(NUM_ELEMS)
            ]
            exp_reverse = list(reversed(src_data[:NUM_ELEMS]))
            cn = tf.next_check(f"{mnemonic} reverse: result")
            tag = f"tc{cn}"
            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_data")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, (t1), {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data_align(sew)
            tf.data_label(f"{tag}_data", format_data_line(src_data, sew))
            tf.data_label(f"{tag}_idx", format_data_line(indices_reverse, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp_reverse, sew))

            # Test 3: gather (splat element 2)
            indices_splat = [2 * elem_bytes] * NUM_ELEMS
            exp_splat = [src_data[2]] * NUM_ELEMS
            cn = tf.next_check(f"{mnemonic} splat: result")
            tag = f"tc{cn}"
            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_data")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, (t1), {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data_align(sew)
            tf.data_label(f"{tag}_data", format_data_line(src_data, sew))
            tf.data_label(f"{tag}_idx", format_data_line(indices_splat, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp_splat, sew))

            # Test 4: masked indexed load (mask=0b0101)
            mask_bits = 0b0101
            vd_init = [0xDD & M(sew)] * NUM_ELEMS
            indices_m = [i * elem_bytes for i in range(NUM_ELEMS)]
            exp_masked = []
            for i in range(NUM_ELEMS):
                if mask_bits & (1 << i):
                    exp_masked.append(src_data[i])
                else:
                    exp_masked.append(vd_init[i])
            cn = tf.next_check(f"{mnemonic} masked: active+inactive")
            tag = f"tc{cn}"
            tf.blank()
            tf.comment(f"Masked indexed load: mask=0b{mask_bits:04b}")
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_vd")
            tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_mask")
            tf.code("vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            tf.code(f"la t1, {tag}_data")
            tf.code(f"{mnemonic} {VREG_DST}, (t1), {VREG_SRC1}, v0.t")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code("la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_vd", format_data_line(vd_init, sew))
            tf.data_label(f"{tag}_data", format_data_line(src_data, sew))
            tf.data_label(f"{tag}_idx", format_data_line(indices_m, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp_masked, sew))

            tf.write(fpath)
            generated.append(str(fpath))

    # --- Indexed stores: vsuxei / vsoxei ---
    for prefix, subdir_name in [("vsuxei", "store"), ("vsoxei", "store")]:
        for sew in SEWS:
            mnemonic = f"{prefix}{sew}.v"
            fname = f"{prefix}{sew}.S"
            fpath = base_dir / "tests" / subdir_name / fname
            tf = TestFile(mnemonic, f"Indexed store {mnemonic}")

            elem_bytes = sew // 8
            src_vals = [10, 20, 30, 40]

            # Test 1: identity store (sequential offsets)
            indices_identity = [i * elem_bytes for i in range(NUM_ELEMS)]
            cn = tf.next_check(f"{mnemonic} identity: result")
            tag = f"tc{cn}"
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            # Load data to store
            tf.code(f"la t1, {tag}_data")
            tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
            # Load indices
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            # Clear result buffer first
            tf.code(f"la t1, result_buf")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, (t1), {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            nbytes = NUM_ELEMS * elem_bytes
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data_align(sew)
            tf.data_label(f"{tag}_data", format_data_line(src_vals, sew))
            tf.data_label(f"{tag}_idx", format_data_line(indices_identity, sew))
            tf.data_label(f"{tag}_exp", format_data_line(src_vals, sew))

            # Test 2: reverse store
            indices_reverse = [
                (NUM_ELEMS - 1 - i) * elem_bytes for i in range(NUM_ELEMS)
            ]
            exp_reverse = list(reversed(src_vals))
            cn = tf.next_check(f"{mnemonic} reverse: result")
            tag = f"tc{cn}"
            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_data")
            tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, result_buf")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} {VREG_DST}, (t1), {VREG_SRC1}")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data_align(sew)
            tf.data_label(f"{tag}_data", format_data_line(src_vals, sew))
            tf.data_label(f"{tag}_idx", format_data_line(indices_reverse, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp_reverse, sew))

            # Test 3: masked indexed store (mask=0b0101)
            mask_bits = 0b0101
            store_data = [0x11, 0x22, 0x33, 0x44]
            store_data = [v & M(sew) for v in store_data]
            buf_init = [0xFF & M(sew)] * NUM_ELEMS
            indices_m = [i * elem_bytes for i in range(NUM_ELEMS)]
            exp_store = []
            for i in range(NUM_ELEMS):
                if mask_bits & (1 << i):
                    exp_store.append(store_data[i])
                else:
                    exp_store.append(buf_init[i])
            cn = tf.next_check(f"{mnemonic} masked: only active stored")
            tag = f"tc{cn}"
            tf.blank()
            tf.comment(f"Masked indexed store: mask=0b{mask_bits:04b}")
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code(f"vle{sew}.v {VREG_DST}, (t1)")
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle{sew}.v {VREG_SRC1}, (t1)")
            tf.code(f"la t1, {tag}_mask")
            tf.code("vlm.v v0, (t1)")
            # Pre-fill result_buf
            tf.code(f"la t1, {tag}_buf_init")
            tf.code(f"vle{sew}.v {VREG_SRC2}, (t1)")
            tf.code("la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_SRC2}, (t1)")
            tf.code("SAVE_CSRS")
            tf.code("la t1, result_buf")
            tf.code(f"{mnemonic} {VREG_DST}, (t1), {VREG_SRC1}, v0.t")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {mask_bits}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_src", format_data_line(store_data, sew))
            tf.data_label(f"{tag}_buf_init", format_data_line(buf_init, sew))
            tf.data_label(f"{tag}_idx", format_data_line(indices_m, sew))
            tf.data_label(f"{tag}_exp", format_data_line(exp_store, sew))

            tf.write(fpath)
            generated.append(str(fpath))

    return generated
