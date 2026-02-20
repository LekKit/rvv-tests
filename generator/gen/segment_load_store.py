from __future__ import annotations

"""Generate tests for segment load/store instructions.

Covers all 252 mnemonics across 9 addressing modes:
  - Unit-stride load:       vlseg{nf}e{eew}.v         (28)
  - Unit-stride store:      vsseg{nf}e{eew}.v         (28)
  - Strided load:           vlsseg{nf}e{eew}.v        (28)
  - Strided store:          vssseg{nf}e{eew}.v        (28)
  - Indexed unordered load: vluxseg{nf}ei{eew}.v      (28)
  - Indexed ordered load:   vloxseg{nf}ei{eew}.v      (28)
  - Indexed unord. store:   vsuxseg{nf}ei{eew}.v      (28)
  - Indexed ordered store:  vsoxseg{nf}ei{eew}.v      (28)
  - Fault-first load:       vlseg{nf}e{eew}ff.v       (28)

nf = 2..8 (7 values), eew = 8/16/32/64 (4 values) → 28 per mode.

Memory layout for a segment of nf fields, vl elements (unit-stride):
    addr[e * nf + f] = field_data[f][e]
After vlseg: register v(8+f) = [field_data[f][0], ..., field_data[f][vl-1]]

Register allocation (LMUL=1):
    v8 .. v(8+nf-1) : destination/source field registers
    v16              : index vector (for indexed variants)

Constraint: EMUL * nf <= 8.  With LMUL=1, EMUL=1, so nf <= 8.  OK.
"""

from pathlib import Path

from ..common import NUM_ELEMS, format_data_line
from ..testfile import TestFile

_NFS = list(range(2, 9))   # nf = 2, 3, ..., 8
_EEWS = [8, 16, 32, 64]


# ------------------------------------------------------------------
# Data generation helpers
# ------------------------------------------------------------------

def _field_data(nf: int) -> list[list[int]]:
    """Test values per field: field_data[f][e] = f * NUM_ELEMS + e + 1."""
    return [[f * NUM_ELEMS + e + 1 for e in range(NUM_ELEMS)] for f in range(nf)]


def _interleaved(nf: int) -> list[int]:
    """Memory layout for unit-stride segment access (interleaved fields)."""
    fields = _field_data(nf)
    mem: list[int] = []
    for e in range(NUM_ELEMS):
        for f in range(nf):
            mem.append(fields[f][e])
    return mem


# ------------------------------------------------------------------
# Unit-stride segment loads: vlseg{nf}e{eew}.v
# ------------------------------------------------------------------

def _gen_seg_loads_unit(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "seg_load"

    for nf in _NFS:
        for eew in _EEWS:
            mnemonic = f"vlseg{nf}e{eew}.v"
            fname = f"vlseg{nf}e{eew}.S"
            tf = TestFile(mnemonic, f"Unit-stride segment load nf={nf} eew={eew}")

            mem = _interleaved(nf)
            fields = _field_data(nf)
            nbytes_per_field = NUM_ELEMS * (eew // 8)

            # Allocate check numbers: one per field + one for CSR
            cns: list[int] = []
            for f in range(nf):
                cns.append(tf.next_check(f"{mnemonic}: field {f} result"))
            cn_csr = tf.next_check(f"{mnemonic}: CSR side-effect")
            tag = f"tc{cns[0]}"

            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{eew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} v8, (t1)")

            # Verify each field register
            for f in range(nf):
                tf.code(f"SET_TEST_NUM {cns[f]}")
                tf.code("la t1, result_buf")
                tf.code(f"vse{eew}.v v{8 + f}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp_f{f}, {nbytes_per_field}")

            tf.code(f"SET_TEST_NUM {cn_csr}")
            tf.code("CHECK_CSRS_UNCHANGED")

            # Data
            tf.data_align(eew)
            tf.data_label(f"{tag}_src", format_data_line(mem, eew))
            for f in range(nf):
                tf.data_label(f"{tag}_exp_f{f}", format_data_line(fields[f], eew))

            fpath = out / fname
            tf.write(fpath)
            generated.append(str(fpath))

    return generated


# ------------------------------------------------------------------
# Unit-stride segment stores: vsseg{nf}e{eew}.v
# ------------------------------------------------------------------

def _gen_seg_stores_unit(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "seg_store"

    for nf in _NFS:
        for eew in _EEWS:
            mnemonic = f"vsseg{nf}e{eew}.v"
            fname = f"vsseg{nf}e{eew}.S"
            tf = TestFile(mnemonic, f"Unit-stride segment store nf={nf} eew={eew}")

            mem_expected = _interleaved(nf)
            fields = _field_data(nf)
            total_bytes = nf * NUM_ELEMS * (eew // 8)

            cn_mem = tf.next_check(f"{mnemonic}: memory result")
            cn_csr = tf.next_check(f"{mnemonic}: CSR side-effect")
            tag = f"tc{cn_mem}"

            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{eew}, m1, tu, mu")

            # Load each field register
            for f in range(nf):
                tf.code(f"la t1, {tag}_f{f}")
                tf.code(f"vle{eew}.v v{8 + f}, (t1)")

            tf.code("SAVE_CSRS")
            tf.code("la t1, result_buf")
            tf.code(f"{mnemonic} v8, (t1)")

            tf.code(f"SET_TEST_NUM {cn_mem}")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {total_bytes}")

            tf.code(f"SET_TEST_NUM {cn_csr}")
            tf.code("CHECK_CSRS_UNCHANGED")

            # Data
            tf.data_align(eew)
            for f in range(nf):
                tf.data_label(f"{tag}_f{f}", format_data_line(fields[f], eew))
            tf.data_label(f"{tag}_exp", format_data_line(mem_expected, eew))

            fpath = out / fname
            tf.write(fpath)
            generated.append(str(fpath))

    return generated


# ------------------------------------------------------------------
# Strided segment loads: vlsseg{nf}e{eew}.v
# ------------------------------------------------------------------

def _gen_seg_loads_strided(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "seg_load"

    for nf in _NFS:
        for eew in _EEWS:
            mnemonic = f"vlsseg{nf}e{eew}.v"
            fname = f"vlsseg{nf}e{eew}.S"
            tf = TestFile(mnemonic, f"Strided segment load nf={nf} eew={eew}")

            # Stride = nf * element_bytes (tightly packed = same as unit-stride)
            seg_bytes = nf * (eew // 8)
            mem = _interleaved(nf)
            fields = _field_data(nf)
            nbytes_per_field = NUM_ELEMS * (eew // 8)

            cns: list[int] = []
            for f in range(nf):
                cns.append(tf.next_check(f"{mnemonic}: field {f} result"))
            cn_csr = tf.next_check(f"{mnemonic}: CSR side-effect")
            tag = f"tc{cns[0]}"

            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{eew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code(f"li t2, {seg_bytes}")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} v8, (t1), t2")

            for f in range(nf):
                tf.code(f"SET_TEST_NUM {cns[f]}")
                tf.code("la t1, result_buf")
                tf.code(f"vse{eew}.v v{8 + f}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp_f{f}, {nbytes_per_field}")

            tf.code(f"SET_TEST_NUM {cn_csr}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(eew)
            tf.data_label(f"{tag}_src", format_data_line(mem, eew))
            for f in range(nf):
                tf.data_label(f"{tag}_exp_f{f}", format_data_line(fields[f], eew))

            fpath = out / fname
            tf.write(fpath)
            generated.append(str(fpath))

    return generated


# ------------------------------------------------------------------
# Strided segment stores: vssseg{nf}e{eew}.v
# ------------------------------------------------------------------

def _gen_seg_stores_strided(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "seg_store"

    for nf in _NFS:
        for eew in _EEWS:
            mnemonic = f"vssseg{nf}e{eew}.v"
            fname = f"vssseg{nf}e{eew}.S"
            tf = TestFile(mnemonic, f"Strided segment store nf={nf} eew={eew}")

            seg_bytes = nf * (eew // 8)
            mem_expected = _interleaved(nf)
            fields = _field_data(nf)
            total_bytes = nf * NUM_ELEMS * (eew // 8)

            cn_mem = tf.next_check(f"{mnemonic}: memory result")
            cn_csr = tf.next_check(f"{mnemonic}: CSR side-effect")
            tag = f"tc{cn_mem}"

            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{eew}, m1, tu, mu")

            for f in range(nf):
                tf.code(f"la t1, {tag}_f{f}")
                tf.code(f"vle{eew}.v v{8 + f}, (t1)")

            tf.code("SAVE_CSRS")
            tf.code("la t1, result_buf")
            tf.code(f"li t2, {seg_bytes}")
            tf.code(f"{mnemonic} v8, (t1), t2")

            tf.code(f"SET_TEST_NUM {cn_mem}")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {total_bytes}")

            tf.code(f"SET_TEST_NUM {cn_csr}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(eew)
            for f in range(nf):
                tf.data_label(f"{tag}_f{f}", format_data_line(fields[f], eew))
            tf.data_label(f"{tag}_exp", format_data_line(mem_expected, eew))

            fpath = out / fname
            tf.write(fpath)
            generated.append(str(fpath))

    return generated


# ------------------------------------------------------------------
# Indexed segment loads: vluxseg / vloxseg
# ------------------------------------------------------------------

def _gen_seg_loads_indexed(base_dir: Path, prefix: str) -> list[str]:
    """Generate indexed segment load tests.

    prefix is 'vluxseg' or 'vloxseg'.
    For each test, SEW (data width) = index EEW for simplicity.
    Indices are sequential: [0, seg_bytes, 2*seg_bytes, 3*seg_bytes].
    """
    generated: list[str] = []
    out = base_dir / "tests" / "seg_load"

    for nf in _NFS:
        for eew in _EEWS:
            mnemonic = f"{prefix}{nf}ei{eew}.v"
            fname = mnemonic.replace(".", "_") + ".S"
            tf = TestFile(mnemonic, f"Indexed segment load nf={nf} eew={eew}")

            seg_bytes = nf * (eew // 8)
            mem = _interleaved(nf)
            fields = _field_data(nf)
            nbytes_per_field = NUM_ELEMS * (eew // 8)

            # Sequential indices (byte offsets)
            indices = [e * seg_bytes for e in range(NUM_ELEMS)]

            cns: list[int] = []
            for f in range(nf):
                cns.append(tf.next_check(f"{mnemonic}: field {f} result"))
            cn_csr = tf.next_check(f"{mnemonic}: CSR side-effect")
            tag = f"tc{cns[0]}"

            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{eew}, m1, tu, mu")

            # Load index vector into v16
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle{eew}.v v16, (t1)")

            tf.code(f"la t1, {tag}_src")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} v8, (t1), v16")

            for f in range(nf):
                tf.code(f"SET_TEST_NUM {cns[f]}")
                tf.code("la t1, result_buf")
                tf.code(f"vse{eew}.v v{8 + f}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp_f{f}, {nbytes_per_field}")

            tf.code(f"SET_TEST_NUM {cn_csr}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(eew)
            tf.data_label(f"{tag}_src", format_data_line(mem, eew))
            tf.data_label(f"{tag}_idx", format_data_line(indices, eew))
            for f in range(nf):
                tf.data_label(f"{tag}_exp_f{f}", format_data_line(fields[f], eew))

            fpath = out / fname
            tf.write(fpath)
            generated.append(str(fpath))

    return generated


# ------------------------------------------------------------------
# Indexed segment stores: vsuxseg / vsoxseg
# ------------------------------------------------------------------

def _gen_seg_stores_indexed(base_dir: Path, prefix: str) -> list[str]:
    """Generate indexed segment store tests.

    prefix is 'vsuxseg' or 'vsoxseg'.
    SEW (data width) = index EEW for simplicity.
    """
    generated: list[str] = []
    out = base_dir / "tests" / "seg_store"

    for nf in _NFS:
        for eew in _EEWS:
            mnemonic = f"{prefix}{nf}ei{eew}.v"
            fname = mnemonic.replace(".", "_") + ".S"
            tf = TestFile(mnemonic, f"Indexed segment store nf={nf} eew={eew}")

            seg_bytes = nf * (eew // 8)
            mem_expected = _interleaved(nf)
            fields = _field_data(nf)
            total_bytes = nf * NUM_ELEMS * (eew // 8)
            indices = [e * seg_bytes for e in range(NUM_ELEMS)]

            cn_mem = tf.next_check(f"{mnemonic}: memory result")
            cn_csr = tf.next_check(f"{mnemonic}: CSR side-effect")
            tag = f"tc{cn_mem}"

            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{eew}, m1, tu, mu")

            for f in range(nf):
                tf.code(f"la t1, {tag}_f{f}")
                tf.code(f"vle{eew}.v v{8 + f}, (t1)")

            # Load index vector into v16
            tf.code(f"la t1, {tag}_idx")
            tf.code(f"vle{eew}.v v16, (t1)")

            tf.code("SAVE_CSRS")
            tf.code("la t1, result_buf")
            tf.code(f"{mnemonic} v8, (t1), v16")

            tf.code(f"SET_TEST_NUM {cn_mem}")
            tf.code(f"CHECK_MEM result_buf, {tag}_exp, {total_bytes}")

            tf.code(f"SET_TEST_NUM {cn_csr}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(eew)
            for f in range(nf):
                tf.data_label(f"{tag}_f{f}", format_data_line(fields[f], eew))
            tf.data_label(f"{tag}_idx", format_data_line(indices, eew))
            tf.data_label(f"{tag}_exp", format_data_line(mem_expected, eew))

            fpath = out / fname
            tf.write(fpath)
            generated.append(str(fpath))

    return generated


# ------------------------------------------------------------------
# Fault-first segment loads: vlseg{nf}e{eew}ff.v
# ------------------------------------------------------------------

def _gen_seg_loads_ff(base_dir: Path) -> list[str]:
    """Generate fault-first segment load tests.

    All data is accessible, so vl should remain unchanged.
    Same verification as unit-stride segment loads.
    """
    generated: list[str] = []
    out = base_dir / "tests" / "seg_load"

    for nf in _NFS:
        for eew in _EEWS:
            mnemonic = f"vlseg{nf}e{eew}ff.v"
            fname = f"vlseg{nf}e{eew}ff.S"
            tf = TestFile(mnemonic, f"Fault-first segment load nf={nf} eew={eew}")

            mem = _interleaved(nf)
            fields = _field_data(nf)
            nbytes_per_field = NUM_ELEMS * (eew // 8)

            cns: list[int] = []
            for f in range(nf):
                cns.append(tf.next_check(f"{mnemonic}: field {f} result"))
            cn_csr = tf.next_check(f"{mnemonic}: CSR side-effect")
            tag = f"tc{cns[0]}"

            tf.code(f"li t0, {NUM_ELEMS}")
            tf.code(f"vsetvli t0, t0, e{eew}, m1, tu, mu")
            tf.code(f"la t1, {tag}_src")
            tf.code("SAVE_CSRS")
            tf.code(f"{mnemonic} v8, (t1)")

            for f in range(nf):
                tf.code(f"SET_TEST_NUM {cns[f]}")
                tf.code("la t1, result_buf")
                tf.code(f"vse{eew}.v v{8 + f}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp_f{f}, {nbytes_per_field}")

            tf.code(f"SET_TEST_NUM {cn_csr}")
            tf.code("CHECK_CSRS_UNCHANGED")

            tf.data_align(eew)
            tf.data_label(f"{tag}_src", format_data_line(mem, eew))
            for f in range(nf):
                tf.data_label(f"{tag}_exp_f{f}", format_data_line(fields[f], eew))

            fpath = out / fname
            tf.write(fpath)
            generated.append(str(fpath))

    return generated


# ------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------

def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []

    # Unit-stride loads/stores
    generated.extend(_gen_seg_loads_unit(base_dir))
    generated.extend(_gen_seg_stores_unit(base_dir))

    # Strided loads/stores
    generated.extend(_gen_seg_loads_strided(base_dir))
    generated.extend(_gen_seg_stores_strided(base_dir))

    # Indexed unordered loads/stores
    generated.extend(_gen_seg_loads_indexed(base_dir, "vluxseg"))
    generated.extend(_gen_seg_loads_indexed(base_dir, "vloxseg"))
    generated.extend(_gen_seg_stores_indexed(base_dir, "vsuxseg"))
    generated.extend(_gen_seg_stores_indexed(base_dir, "vsoxseg"))

    # Fault-first loads
    generated.extend(_gen_seg_loads_ff(base_dir))

    return generated
