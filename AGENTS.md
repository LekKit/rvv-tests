# AI-Generated Code Declaration

This project was entirely generated through AI-assisted development (vibecoding).

## Tools used

- **[OpenCode](https://github.com/anomalyco/opencode)** with Anthropic Claude Opus 4.6 — all code generation, architecture design, debugging, and test vector computation
- **Google Gemini 3.1 Pro** (chat, not agentic) — architectural review and test strategy suggestions: register overlap hazard testing, VLEN-agnostic tail checking via `vs1r.v`/`vlenb`, extreme arithmetic corner cases, LMUL>1 boundary crossing, mask/tail policy interaction, vstart platform limitations diagnosis, vl=0 no-op coverage, fault-only-first strategy
- **Human operator** — project direction, review, and deployment coordination

## What was AI-generated

Everything:
- Project architecture and design decisions
- Python test generator framework (`generator/`)
- Assembly test macros and infrastructure (`include/`)
- All compute functions for expected value calculation
- Build system (`Makefile`)
- Test runner scripts (`scripts/`)
- Documentation

## Project structure

```
rvv-tests/
├── generate_tests.py              # Entry point: python3 generate_tests.py -v
├── Makefile                       # Cross-compilation build system
├── include/
│   ├── riscv_test.h               # _start entry, _mem_compare, stack setup
│   └── test_macros.h              # CHECK_MEM, SAVE_CSRS, FAIL_IF_*, witness macros, SYS_* syscall wrappers
├── scripts/
│   └── run_tests.sh               # Test runner (local + SSH modes)
├── generator/                     # Python test generator framework
│   ├── main.py                    # Driver: imports generators, cleans tests/, writes MANIFEST
│   ├── common.py                  # Shared constants (SEWS, NUM_ELEMS, M/S/U, IEEE754)
│   ├── testfile.py                # TestFile class (assembly file builder)
│   ├── emit.py                    # Reusable emit helpers (binop, widening, compare, ...)
│   ├── vectors.py                 # Test vector definitions per instruction family
│   ├── compute/                   # Expected-value computation
│   │   ├── integer.py             # Integer arithmetic, shifts, compare, div, MAC
│   │   ├── fixed_point.py         # Saturating, averaging, scaling, narrowing clip
│   │   └── floating_point.py      # FP arith, conversions, frsqrt7/frec7 spec tables
│   └── gen/                       # Per-family test generators (one file per family)
│       ├── config.py              # vsetvli / vsetivli / vsetvl
│       ├── load_store.py          # Unit-stride, strided, whole-register, fault-first
│       ├── indexed_load_store.py  # vluxei / vloxei / vsuxei / vsoxei
│       ├── segment_load_store.py  # All segment load/store variants (252 tests)
│       ├── int_arith.py           # add/sub/logical/shift/minmax/mul/div
│       ├── int_cmp.py             # Integer compares (vmseq..vmsgt)
│       ├── int_widening.py        # Widening arith + narrowing shifts
│       ├── int_macc.py            # Integer + widening multiply-accumulate
│       ├── int_extension.py       # vzext / vsext
│       ├── int_adc.py             # Add/subtract with carry/borrow
│       ├── fixed_point.py         # All fixed-point operations
│       ├── float_arith.py         # FP arith/cmp/convert/misc/rsqrt7/rec7
│       ├── float_widening.py      # Widening/narrowing FP arith + conversions
│       ├── reduction.py           # Integer + FP reductions
│       ├── mask.py                # Mask logical + population/find/set/iota/vid
│       ├── permutation.py         # Slides, gathers, compress, moves, merges
│       └── edge_cases.py          # vl=0, tail undisturbed, LMUL>1, fault-only-first, vill trap, fractional LMUL, tail/mask agnostic, stride, scatter, vxsat sticky, narrowing/widening
├── tests/                         # Generated .S files (653 total) — fully regenerated on each run
│   ├── MANIFEST                   # List of all test file paths
│   ├── config/                    # vsetvli (covers vsetvli/vsetivli/vsetvl)
│   ├── load/                      # vle*, vlse*, vlm, vle*ff, vl*re*
│   ├── store/                     # vse*, vsse*, vs*r
│   ├── seg_load/                  # vlseg*, vlsseg*, vluxseg*, vloxseg*, vlseg*ff
│   ├── seg_store/                 # vsseg*, vssseg*, vsuxseg*, vsoxseg*
│   ├── int_arith/                 # vadd, vsub, vrsub
│   ├── int_logical/               # vand, vor, vxor
│   ├── int_shift/                 # vsll, vsrl, vsra
│   ├── int_minmax/                # vmin[u], vmax[u]
│   ├── int_mul/                   # vmul, vmulh[su][u]
│   ├── int_div/                   # vdiv[u], vrem[u]
│   ├── int_cmp/                   # vmseq..vmsgt
│   ├── int_widening/              # vwadd..vwmulsu, vnsrl, vnsra
│   ├── int_macc/                  # vmacc..vwmaccus
│   ├── int_extension/             # vzext, vsext
│   ├── int_adc/                   # vadc, vsbc, vmadc, vmsbc
│   ├── fixed_point/               # vsadd..vnclip
│   ├── float_arith/               # vfadd..vfrdiv
│   ├── float_cmp/                 # vmfeq..vmfge
│   ├── float_convert/             # vfcvt.*
│   ├── float_minmax/              # vfmin, vfmax
│   ├── float_misc/                # vfsqrt, vfclass, vfrsqrt7, vfrec7
│   ├── float_muladd/              # vfmacc..vfnmsub
│   ├── float_sgnj/                # vfsgnj[n][x]
│   ├── float_widening/            # vfwadd..vfwcvt*
│   ├── float_narrowing/           # vfncvt* (including rod)
│   ├── reduction/                 # vredsum..vfwredusum
│   ├── mask/                      # vmand..vid
│   ├── permutation/               # vslide..vmv*r
│   └── edge_cases/                # vl=0, tail, LMUL>1, vle32ff fault, vill trap, fractional LMUL, tail/mask agnostic, stride, scatter, vxsat sticky, narrowing/widening
└── fuzzer/                        # Planned AFL++ fuzzing (not yet implemented)
```

## Internal design details

### Test format

Each `.S` file is a **standalone static binary** — no libc, no dynamic linking. The entry point `_start` (from `riscv_test.h`) sets up a small stack and jumps to `_test_start`. Tests exit via Linux syscall 93: exit code 0 means all checks passed, any nonzero code identifies the first check that failed. Check meanings are documented in comments at the top of each generated `.S` file.

Edge case tests also use Linux syscalls 222 (mmap) and 215 (munmap) to create page-boundary conditions for fault-only-first load testing, and syscall 220 (clone) + 260 (wait4) to fork child processes for illegal-instruction trap verification.

### Register allocation convention

All tests use a fixed register allocation to keep things predictable:

| Register | Role |
|----------|------|
| `v8` | Primary destination (VREG_DST) |
| `v16` | vs2 / first vector source (VREG_SRC2) |
| `v20` | vs1 / second vector source (VREG_SRC1) |
| `v24` | Witness register (should not be modified) |
| `v0` | Mask register (when needed) |
| `a0` | Scalar operand for VX forms |
| `fa0` | FP scalar operand for VF forms |
| `s0`–`s5` | Saved CSR snapshots (vstart, vxsat, vxrm, vl, vtype, fcsr) |
| `s6`–`s10` | Free for test-specific use (e.g. mmap base, load pointer, saved vl) |
| `s11` | Current test/check number (used as exit code on failure) |

### VLEN-agnostic design

Tests use `vsetvli` with `vl=4` and `LMUL=1`, so they work on any hardware with `VLEN >= 128`. The fixed small element count (NUM_ELEMS=4) keeps test data manageable and avoids VLEN-dependent behavior.

Full-VLMAX tail tests are the exception: they use `vsetvli t0, x0, ...` (request VLMAX), dump with `vs1r.v`/`vs2r.v` (stores exactly VLEN/8 or 2×VLEN/8 bytes regardless of vtype/vl), and compute the tail region size at runtime via `csrr vlenb`.

### Side-effect checking

After executing the instruction under test, each test verifies:

1. **vstart == 0** — all vector instructions must reset vstart
2. **vxsat unchanged** — unless the instruction is a fixed-point saturating operation
3. **vxrm unchanged** — always (no instruction writes vxrm as a side effect)
4. **vl unchanged** — unless the instruction is a configuration instruction
5. **vtype unchanged** — unless the instruction is a configuration instruction
6. **fcsr unchanged** — unless the instruction is a floating-point operation
7. **Witness register unchanged** — v24 is pre-loaded with a known pattern and verified after the instruction, catching writes to wrong destination registers

Three macro variants handle different instruction classes:
- `CHECK_CSRS_UNCHANGED` — full check (integer/load/store/mask/permutation)
- `CHECK_CSRS_UNCHANGED_FP` — skips fcsr (FP ops may set fflags)
- `CHECK_CSRS_UNCHANGED_FIXEDPOINT` — skips vxsat (saturating ops may set it)

### Platform limitations

**Writing vstart from userspace causes SIGILL** on Linux targets where the kernel does not support trap-and-emulate for vstart writes. This is legal per the RVV spec: hardware that never interrupts vector instructions mid-execution need not support arbitrary vstart values, and the kernel may trap writes to vstart as a high-cost operation. The `vstart_nonzero` test is therefore commented out (not deleted) in `edge_cases.py`; it is intended for a future bare-metal port.

**Userspace SIGILL handler PC-advance is unreliable**: The standard approach of using `rt_sigaction` with `SA_SIGINFO` to catch SIGILL and advance the saved PC (at ucontext offset 168 on RV64) does not work reliably across all platforms/emulators. The `vill_trap` test instead uses `clone`/`wait4` to fork a child process that executes the faulting instruction; the parent checks the child was killed by SIGILL (WTERMSIG == 4).

**Vector state (vtype/vl) is not preserved across fork**: The Linux kernel's lazy vector context management may not preserve `vtype` and `vl` CSRs across `clone()`/`fork()`. Child processes must set their own vector configuration rather than relying on inheriting the parent's vtype. Vector register contents may also not survive `clone`+`wait4` syscall sequences in the parent.

### Code generation architecture

The generator follows a pipeline:

1. **`vectors.py`** defines test input vectors per instruction family (e.g., `binop_vv()` yields named tuples of source operands for each SEW)
2. **`compute/*.py`** contains pure-Python functions that compute the expected result for each instruction (e.g., `add(a, b, sew)` returns the SEW-truncated sum)
3. **`gen/*.py`** files iterate over SEW widths and test vectors, calling compute functions to get expected values, then calling emit helpers to generate the assembly
4. **`emit.py`** provides reusable templates (`emit_binop_vv`, `emit_widening_vv`, `emit_compare_vv`, etc.) that handle the boilerplate: setup vtype, load sources, save CSRs, execute instruction, store result, compare, check CSRs, check witness
5. **`testfile.py`** (`TestFile` class) accumulates `.text` and `.data` sections and writes the final `.S` file with proper includes and check number comments

On each run, `main.py` removes the entire `tests/` directory before regenerating (safety: only proceeds if the directory is named `tests`), ensuring no stale files from removed tests accumulate.

### Chicken-and-egg resolution

Load/store tests are verified using scalar instructions (`lbu`/`lhu`/`lwu`/`ld`) to check memory contents — these are in the GC extension which is assumed correct. All arithmetic tests assume loads/stores work (since those are tested first).

### Floating-point specifics

- **frsqrt7.v / frec7.v**: These are approximate instructions using spec-defined lookup tables (RVV 1.0 Tables 16 and 17). The compute functions implement the exact spec algorithm with the full lookup tables, not a math library approximation.
- **vfncvt.rod.f.f.w**: Uses round-to-odd semantics. Test vectors use exact values where all rounding modes agree.
- **FP compute uses Python f64**: f32 results are obtained via `struct.pack/unpack` round-trip. Test vectors are chosen to avoid double-rounding issues.
- **`math.fma`**: Only available in Python 3.13+; the code has a `hasattr` fallback that uses `a * b + c` for older versions.

### Buffer sizes

- `result_buf`: 256 bytes — stores instruction output for memory comparison
- `witness_buf`: 256 bytes — stores witness register for verification
- Maximum usage: `nf=8, eew=64, 4 elements` = 8 * 4 * 8 = 256 bytes (exactly at capacity for segment load/store tests)
- Full-VLMAX and widening edge case tests use their own larger buffers (512–1024 bytes) in `.data`, sized to support up to VLEN=4096
