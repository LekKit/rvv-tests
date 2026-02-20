from __future__ import annotations

"""Generate tests for integer extension instructions:
vzext.vf2, vzext.vf4, vzext.vf8,
vsext.vf2, vsext.vf4, vsext.vf8.
"""

from pathlib import Path

from ..common import SEWS, NUM_ELEMS, M, S, U, format_data_line, format_hex
from ..testfile import TestFile
from ..emit import VREG_DST, VREG_SRC2, VREG_WITNESS
from ..compute.integer import zext2, sext2, zext4, sext4, zext8, sext8

# (mnemonic, dest_sew_list, factor, compute_fn)
_EXT_OPS: list[tuple[str, list[int], int, object]] = [
    ("vzext.vf2", [16, 32, 64], 2, zext2),
    ("vzext.vf4", [32, 64],     4, zext4),
    ("vzext.vf8", [64],         8, zext8),
    ("vsext.vf2", [16, 32, 64], 2, sext2),
    ("vsext.vf4", [32, 64],     4, sext4),
    ("vsext.vf8", [64],         8, sext8),
]


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "int_extension"

    for mnemonic, dest_sews, factor, cfn in _EXT_OPS:
        fname = mnemonic.replace(".", "_") + ".S"
        tf = TestFile(mnemonic, f"Integer extension {mnemonic}")

        for dsew in dest_sews:
            ssew = dsew // factor  # source element width
            sm = M(ssew)
            # Test vectors at source width
            test_cases = [
                ("basic", [1, 2, 3, 4]),
                ("max",   [sm, sm - 1, sm >> 1, 0]),
                ("high",  [sm, (1 << (ssew - 1)), (1 << (ssew - 1)) + 1, 1]),
            ]
            for name, src in test_cases:
                exp = [cfn(v, dsew) for v in src]
                nbytes = NUM_ELEMS * (dsew // 8)
                cn_res = tf.next_check(
                    f"{mnemonic} e{ssew}→e{dsew} {name}: result")
                cn_csr = tf.next_check(
                    f"{mnemonic} e{ssew}→e{dsew} {name}: CSR side-effect")
                tag = f"tc{cn_res}"

                tf.blank()
                tf.comment(
                    f"Test {cn_res}-{cn_csr}: {mnemonic} "
                    f"e{ssew}→e{dsew} {name}")
                # Load source at source SEW
                tf.code(f"li t0, {NUM_ELEMS}")
                tf.code(f"vsetvli t0, t0, e{ssew}, m1, tu, mu")
                tf.code(f"la t1, {tag}_src")
                tf.code(f"vle{ssew}.v {VREG_SRC2}, (t1)")
                # Execute with dest SEW set (LMUL adjusted)
                lmul_d = max(1, factor)  # dest LMUL = factor
                tf.code(f"vsetvli t0, t0, e{dsew}, m{lmul_d}, tu, mu")
                tf.code("SAVE_CSRS")
                tf.code(f"{mnemonic} {VREG_DST}, {VREG_SRC2}")
                # Store result at dest SEW
                tf.code(f"SET_TEST_NUM {cn_res}")
                tf.code(f"la t1, result_buf")
                tf.code(f"vse{dsew}.v {VREG_DST}, (t1)")
                tf.code(f"CHECK_MEM result_buf, {tag}_exp, {nbytes}")
                tf.code(f"SET_TEST_NUM {cn_csr}")
                tf.code("CHECK_VSTART_ZERO")

                tf.data_align(dsew)
                tf.data_label(f"{tag}_src", format_data_line(src, ssew))
                tf.data_label(f"{tag}_exp", format_data_line(exp, dsew))

        fpath = out / fname
        tf.write(fpath)
        generated.append(str(fpath))

    return generated
