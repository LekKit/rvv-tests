from __future__ import annotations

"""TestFile — builder for a single self-contained ``.S`` test file."""

from pathlib import Path


class TestFile:
    """Accumulates assembly code + data and writes a complete ``.S`` file.

    Typical usage::

        tf = TestFile("vadd.vv", "Vector-vector integer add")
        # … add test cases via emit helpers …
        tf.write("tests/int_arith/vadd_vv.S")
    """

    def __init__(self, mnemonic: str, description: str) -> None:
        self.mnemonic = mnemonic
        self.description = description
        self._code: list[str] = []
        self._data: list[str] = []
        self._check_num: int = 0
        self._check_comments: list[str] = []

    # ------------------------------------------------------------------
    # Check-number allocator
    # ------------------------------------------------------------------

    @property
    def check_count(self) -> int:
        return self._check_num

    def next_check(self, comment: str = "") -> int:
        """Allocate the next sequential check number (1-based)."""
        self._check_num += 1
        if comment:
            self._check_comments.append(f" *   {self._check_num:3d} = {comment}")
        return self._check_num

    # ------------------------------------------------------------------
    # Code helpers
    # ------------------------------------------------------------------

    def code(self, line: str, indent: int = 1) -> None:
        self._code.append("    " * indent + line)

    def blank(self) -> None:
        self._code.append("")

    def comment(self, text: str, indent: int = 1) -> None:
        self._code.append("    " * indent + f"/* {text} */")

    def raw(self, text: str) -> None:
        """Append raw text (no indentation)."""
        self._code.append(text)

    # ------------------------------------------------------------------
    # Data helpers
    # ------------------------------------------------------------------

    def data(self, line: str) -> None:
        self._data.append(line)

    def data_align(self, sew: int) -> None:
        align = max(2, sew // 8)       # at least 4-byte aligned
        self._data.append(f".align {align.bit_length() - 1}")

    def data_label(self, label: str, content: str) -> None:
        self._data.append(f"{label}:")
        self._data.append(content)

    # ------------------------------------------------------------------
    # Write to disk
    # ------------------------------------------------------------------

    def write(self, filepath: str | Path) -> int:
        """Write the test file and return the total check count."""
        p = Path(filepath)
        p.parent.mkdir(parents=True, exist_ok=True)

        lines: list[str] = []

        # ── header ──
        lines.append(f"/* Auto-generated test for {self.mnemonic}")
        lines.append(f" * {self.description}")
        lines.append(" *")
        lines.append(" * Exit code 0 = PASS")
        lines.append(" * Exit code N = check N failed:")
        lines.extend(self._check_comments)
        lines.append(" */")
        lines.append('#include "riscv_test.h"')
        lines.append('#include "test_macros.h"')
        lines.append("")

        # ── code section ──
        lines.extend(self._code)
        lines.append("")
        lines.append("    PASS_TEST")

        # ── data section ──
        lines.append("")
        lines.append(".data")
        lines.extend(self._data)
        lines.append("")
        lines.append(".align 4")
        lines.append("result_buf:  .space 256")
        lines.append("witness_buf: .space 256")
        lines.append("")

        p.write_text("\n".join(lines) + "\n")
        return self._check_num
