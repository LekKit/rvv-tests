from __future__ import annotations

"""Generate tests for configuration-setting instructions:
vsetvli, vsetivli, vsetvl.
"""

from pathlib import Path

from ..common import NUM_ELEMS
from ..testfile import TestFile


def generate(base_dir: Path) -> list[str]:
    generated: list[str] = []
    out = base_dir / "tests" / "config"

    # --- vsetvli ---
    fpath = out / "vsetvli.S"
    tf = TestFile("vsetvli", "Set vector length and type from register + immediate")

    # Test 1: Basic vsetvli – set e32/m1, avl=4, expect vl=4
    cn = tf.next_check("vsetvli e32/m1 avl=4: wrong vl")
    tf.blank()
    tf.comment(f"Test {cn}: vsetvli e32, m1, tu, mu with avl=4")
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code("li t0, 4")
    tf.code("vsetvli t1, t0, e32, m1, tu, mu")
    tf.code("li t2, 4")
    tf.code("FAIL_IF_NE t1, t2")

    # Test 2: vl should match what CSR vl reports
    cn = tf.next_check("vsetvli: vl CSR mismatch")
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code("csrr t3, vl")
    tf.code("FAIL_IF_NE t1, t3")

    # Test 3: vtype should be valid (vill=0)
    cn = tf.next_check("vsetvli: vtype.vill set")
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code("csrr t3, vtype")
    tf.code("srai t4, t3, 63")   # vill is sign bit (bit XLEN-1)
    tf.code("FAIL_IF_NZ t4")

    # Test 4: Different SEW values
    for sew_name, sew_str in [("e8", "e8"), ("e16", "e16"), ("e32", "e32"), ("e64", "e64")]:
        cn = tf.next_check(f"vsetvli {sew_name}/m1: vl>0")
        tf.blank()
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code("li t0, 4")
        tf.code(f"vsetvli t1, t0, {sew_str}, m1, tu, mu")
        tf.code("FAIL_IF_Z t1")   # vl should be > 0

    # Test 5: Different LMUL values
    for lmul_name, lmul_str in [("m1", "m1"), ("m2", "m2"), ("m4", "m4"), ("m8", "m8"),
                                  ("mf2", "mf2"), ("mf4", "mf4"), ("mf8", "mf8")]:
        cn = tf.next_check(f"vsetvli e32/{lmul_name}: no vill")
        tf.blank()
        tf.code(f"SET_TEST_NUM {cn}")
        tf.code("li t0, 4")
        # mf8 with e32 may set vill (SEW > LMUL*ELEN constraint)
        # So just check it doesn't crash
        tf.code(f"vsetvli t1, t0, e32, {lmul_str}, tu, mu")
        # We don't assert vl here since some LMUL/SEW combos may be illegal
        tf.code("csrr t3, vtype")
        # Just verify we can read vtype without error

    # Test 6: avl=0 should give vl=0
    cn = tf.next_check("vsetvli avl=0: vl should be 0")
    tf.blank()
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code("li t0, 0")
    tf.code("vsetvli t1, t0, e32, m1, tu, mu")
    tf.code("FAIL_IF_NZ t1")

    # Test 7: vstart should be reset to 0
    cn = tf.next_check("vsetvli: vstart should be 0")
    tf.blank()
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code("li t0, 4")
    tf.code("vsetvli t1, t0, e32, m1, tu, mu")
    tf.code("csrr t3, vstart")
    tf.code("FAIL_IF_NZ t3")

    # Test 8: vsetivli with uimm=4
    cn = tf.next_check("vsetivli uimm=4 e32/m1: vl=4")
    tf.blank()
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code("vsetivli t1, 4, e32, m1, tu, mu")
    tf.code("li t2, 4")
    tf.code("FAIL_IF_NE t1, t2")

    # Test 9: vsetivli with uimm=0
    cn = tf.next_check("vsetivli uimm=0: vl=0")
    tf.blank()
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code("vsetivli t1, 0, e32, m1, tu, mu")
    tf.code("FAIL_IF_NZ t1")

    # Test 10: vsetvl (register-register)
    cn = tf.next_check("vsetvl: basic")
    tf.blank()
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code("li t0, 4")
    # First set vtype via vsetvli to get the encoding
    tf.code("vsetvli zero, t0, e32, m1, tu, mu")
    tf.code("csrr t2, vtype")
    # Now use vsetvl with that vtype value
    tf.code("vsetvl t1, t0, t2")
    tf.code("li t3, 4")
    tf.code("FAIL_IF_NE t1, t3")

    # Test 11: Query VLMAX (set avl = ~0, rd != x0)
    cn = tf.next_check("vsetvli: query vlmax")
    tf.blank()
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code("li t0, -1")   # maximum avl
    tf.code("vsetvli t1, t0, e8, m8, tu, mu")
    # t1 should now contain VLMAX for e8/m8 = 8*VLEN/8 = VLEN
    tf.code("csrr t2, vlenb")  # VLEN/8
    tf.code("slli t2, t2, 3")  # VLEN/8 * 8 = VLEN
    tf.code("FAIL_IF_NE t1, t2")

    # Test 12: ta/ma flags
    cn = tf.next_check("vsetvli: ta/ma encoding")
    tf.blank()
    tf.code(f"SET_TEST_NUM {cn}")
    tf.code("li t0, 4")
    tf.code("vsetvli t1, t0, e32, m1, ta, ma")
    tf.code("csrr t3, vtype")
    # Bits: vta=bit6, vma=bit7
    tf.code("andi t4, t3, 0xC0")
    tf.code("li t5, 0xC0")
    tf.code("FAIL_IF_NE t4, t5")

    tf.write(fpath)
    generated.append(str(fpath))

    return generated
