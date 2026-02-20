from __future__ import annotations

"""Generate tests for mask instructions:
vmand, vmnand, vmandn, vmxor, vmor, vmnor, vmorn, vmxnor,
vcpop.m, vfirst.m, vmsbf.m, vmsif.m, vmsof.m,
viota.m, vid.v.
"""

from pathlib import Path
from typing import Callable

from ..common import SEWS, NUM_ELEMS, M, U, format_data_line
from ..testfile import TestFile
from ..emit import emit_mask_logical, VREG_DST, VREG_SRC2, VREG_SRC1
from ..vectors import mask_logical
from ..compute.integer import (
    vmand, vmnand, vmandn, vmxor, vmor, vmnor, vmorn, vmxnor,
)

_MASK_LOGIC: list[tuple[str, Callable]] = [
    ("vmand.mm",  vmand),
    ("vmnand.mm", vmnand),
    ("vmandn.mm", vmandn),
    ("vmxor.mm",  vmxor),
    ("vmor.mm",   vmor),
    ("vmnor.mm",  vmnor),
    ("vmorn.mm",  vmorn),
    ("vmxnor.mm", vmxnor),
]


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "mask"

    # --- Mask logical ---
    for mnemonic, cfn in _MASK_LOGIC:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Mask logical {mnemonic}")
        for name, m2, m1 in mask_logical():
            exp = cfn(m2, m1) & 0xF  # only low 4 bits matter for 4 elements
            emit_mask_logical(tf, mnemonic, name, m2, m1, exp)
        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    # --- vcpop.m ---
    mnemonic = "vcpop.m"
    fpath = out / "vcpop_m.S"
    tf = TestFile(mnemonic, "Mask population count")
    test_cases = [
        ("none",  0b0000, 0),
        ("one",   0b0001, 1),
        ("two",   0b0101, 2),
        ("all",   0b1111, 4),
        ("three", 0b1110, 3),
    ]
    for name, mask_val, expected_count in test_cases:
        cn = tf.next_check(f"vcpop.m {name}: wrong count")
        tag = f"tc{cn}"
        tf.blank()
        tf.comment(f"Test {cn}: vcpop.m {name}")
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e8, m1, tu, mu")
        tf.code(f"la t1, {tag}_mask")
        tf.code(f"vlm.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vcpop.m t2, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"li t3, {expected_count}")
        tf.code("FAIL_IF_NE t2, t3")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data(".align 1")
        tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vfirst.m ---
    mnemonic = "vfirst.m"
    fpath = out / "vfirst_m.S"
    tf = TestFile(mnemonic, "Find first set mask bit")
    test_cases = [
        ("none",  0b0000, -1),
        ("bit0",  0b0001, 0),
        ("bit1",  0b0010, 1),
        ("bit2",  0b0100, 2),
        ("bit3",  0b1000, 3),
        ("multi", 0b1010, 1),
    ]
    for name, mask_val, expected_idx in test_cases:
        cn = tf.next_check(f"vfirst.m {name}: wrong index")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e8, m1, tu, mu")
        tf.code(f"la t1, {tag}_mask")
        tf.code(f"vlm.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vfirst.m t2, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"li t3, {expected_idx}")
        tf.code("FAIL_IF_NE t2, t3")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data(".align 1")
        tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vid.v ---
    mnemonic = "vid.v"
    fpath = out / "vid_v.S"
    tf = TestFile(mnemonic, "Vector element index")
    for sew in SEWS:
        expected = list(range(NUM_ELEMS))
        nbytes = NUM_ELEMS * (sew // 8)
        cn = tf.next_check(f"vid.v e{sew}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
        tf.code("SAVE_CSRS")
        tf.code(f"vid.v {VREG_DST}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
        tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data_align(sew)
        tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- viota.m ---
    mnemonic = "viota.m"
    fpath = out / "viota_m.S"
    tf = TestFile(mnemonic, "Vector iota (prefix sum of mask)")
    for sew in SEWS:
        test_cases_iota = [
            ("all",   0b1111, [0, 1, 2, 3]),
            ("none",  0b0000, [0, 0, 0, 0]),
            ("alt",   0b1010, [0, 0, 1, 1]),
            ("first", 0b0001, [0, 1, 1, 1]),
        ]
        for name, mask_val, expected in test_cases_iota:
            nbytes = NUM_ELEMS * (sew // 8)
            cn = tf.next_check(f"viota.m e{sew} {name}: result")
            tag = f"tc{cn}"
            tf.blank()
            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{sew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_mask")
            tf.code(f"vlm.v v0, (t1)")
            tf.code("SAVE_CSRS")
            # viota.m vd, vs2  — reads mask from vs2
            tf.code(f"viota.m {VREG_DST}, v0")
            tf.code(f"SET_TEST_NUM {cn}")
            tf.code(f"la t1, result_buf")
            tf.code(f"vse{sew}.v {VREG_DST}, (t1)")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
            tf.code("CHECK_CSRS_UNCHANGED")
            tf.data(".align 1")
            tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
            tf.data_align(sew)
            tf.data_label(f"{tag}_exp", format_data_line(expected, sew))
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmsbf.m (set-before-first) ---
    mnemonic = "vmsbf.m"
    fpath = out / "vmsbf_m.S"
    tf = TestFile(mnemonic, "Set before first set mask bit")
    # vmsbf: result[i]=1 for all positions before the first set bit in vs2
    test_cases_sbf = [
        ("none",   0b0000, 0b1111),  # no set bit → all before
        ("bit0",   0b0001, 0b0000),  # first at 0 → nothing before
        ("bit1",   0b0010, 0b0001),  # first at 1 → bit 0
        ("bit2",   0b0100, 0b0011),  # first at 2 → bits 0,1
        ("bit3",   0b1000, 0b0111),  # first at 3 → bits 0,1,2
        ("multi",  0b1010, 0b0001),  # first at 1 → bit 0
    ]
    for name, mask_val, exp_mask in test_cases_sbf:
        cn = tf.next_check(f"vmsbf.m {name}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e8, m1, tu, mu")
        tf.code(f"la t1, {tag}_mask")
        tf.code(f"vlm.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vmsbf.m {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vsm.v {VREG_DST}, (t1)")
        tf.code(f"lbu t2, 0(t1)")
        tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
        tf.code(f"li t3, {exp_mask}")
        tf.code("FAIL_IF_NE t2, t3")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data(".align 1")
        tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmsif.m (set-including-first) ---
    mnemonic = "vmsif.m"
    fpath = out / "vmsif_m.S"
    tf = TestFile(mnemonic, "Set including first set mask bit")
    test_cases_sif = [
        ("none",   0b0000, 0b1111),  # no set bit → all set
        ("bit0",   0b0001, 0b0001),  # first at 0 → just bit 0
        ("bit1",   0b0010, 0b0011),  # first at 1 → bits 0,1
        ("bit2",   0b0100, 0b0111),  # first at 2 → bits 0,1,2
        ("bit3",   0b1000, 0b1111),  # first at 3 → all
        ("multi",  0b1010, 0b0011),  # first at 1 → bits 0,1
    ]
    for name, mask_val, exp_mask in test_cases_sif:
        cn = tf.next_check(f"vmsif.m {name}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e8, m1, tu, mu")
        tf.code(f"la t1, {tag}_mask")
        tf.code(f"vlm.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vmsif.m {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vsm.v {VREG_DST}, (t1)")
        tf.code(f"lbu t2, 0(t1)")
        tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
        tf.code(f"li t3, {exp_mask}")
        tf.code("FAIL_IF_NE t2, t3")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data(".align 1")
        tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
    tf.write(fpath)
    generated.append(str(fpath))

    # --- vmsof.m (set-only-first) ---
    mnemonic = "vmsof.m"
    fpath = out / "vmsof_m.S"
    tf = TestFile(mnemonic, "Set only first set mask bit")
    test_cases_sof = [
        ("none",   0b0000, 0b0000),  # no set bit → none
        ("bit0",   0b0001, 0b0001),  # first at 0
        ("bit1",   0b0010, 0b0010),  # first at 1
        ("bit2",   0b0100, 0b0100),  # first at 2
        ("bit3",   0b1000, 0b1000),  # first at 3
        ("multi",  0b1010, 0b0010),  # first at 1
        ("all",    0b1111, 0b0001),  # first at 0
    ]
    for name, mask_val, exp_mask in test_cases_sof:
        cn = tf.next_check(f"vmsof.m {name}: result")
        tag = f"tc{cn}"
        tf.blank()
        tf.code(f"li t0, {NUM_ELEMS}")
        tf.code(f"vsetvli t0, t0, e8, m1, tu, mu")
        tf.code(f"la t1, {tag}_mask")
        tf.code(f"vlm.v {VREG_SRC2}, (t1)")
        tf.code("SAVE_CSRS")
        tf.code(f"vmsof.m {VREG_DST}, {VREG_SRC2}")
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code(f"la t1, result_buf")
        tf.code(f"vsm.v {VREG_DST}, (t1)")
        tf.code(f"lbu t2, 0(t1)")
        tf.code(f"andi t2, t2, {(1 << NUM_ELEMS) - 1}")
        tf.code(f"li t3, {exp_mask}")
        tf.code("FAIL_IF_NE t2, t3")
        tf.code("CHECK_CSRS_UNCHANGED")
        tf.data(".align 1")
        tf.data_label(f"{tag}_mask", f"    .byte {mask_val}")
    tf.write(fpath)
    generated.append(str(fpath))

    return generated
