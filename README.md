# rvv-tests

Comprehensive test suite for the RISC-V Vector Extension 1.0 (RVV 1.0).

These tests are 100% vibecoded.

## Overview

653 standalone assembly tests covering all RVV 1.0 instruction mnemonics. Each test is a static binary (no libc) that exits with code 0 on success or a nonzero check number identifying the first failed verification.

Tests verify:
- Correct instruction output for known inputs across all supported SEW widths
- No side effects on unrelated V-extension state (vstart, vxsat, vxrm, vl, vtype, fcsr)
- No corruption of bystander vector registers
- Tail-undisturbed policy (elements >= vl preserved) up to full VLEN/8 bytes
- LMUL>1 register boundary crossing correctness
- vl=0 no-op behavior (integer, FP, load, store)
- Fault-only-first load semantics via mmap/munmap page boundary

### Instruction families covered

| Family | Tests | Mnemonics |
|--------|-------|-----------|
| Configuration | 1 | vsetvli, vsetivli, vsetvl |
| Unit-stride load/store | 42 | vle/vse/vlse/vsse/vlm/vl*re/vs*r/vle*ff |
| Indexed load/store | 16 | vluxei/vloxei/vsuxei/vsoxei |
| Segment load/store | 252 | vlseg/vsseg/vlsseg/vssseg/vluxseg/vloxseg/vsuxseg/vsoxseg/vlseg*ff |
| Integer arithmetic | 49 | vadd/vsub/vrsub/vand/vor/vxor/shifts/minmax/mul/div |
| Integer compare | 20 | vmseq/vmsne/vmsltu/vmslt/vmsleu/vmsle/vmsgtu/vmsgt |
| Integer widening/narrowing | 28 | vwadd/vwsub/vwmul/vnsrl/vnsra |
| Integer MAC | 15 | vmacc/vmadd/vnmsac/vnmsub/vwmacc* |
| Integer extension | 6 | vzext/vsext .vf2/.vf4/.vf8 |
| Integer add-with-carry | 15 | vadc/vsbc/vmadc/vmsbc |
| Fixed-point | 32 | vsadd/vssub/vaadd/vasub/vsmul/vssrl/vssra/vnclipu/vnclip |
| FP arithmetic | 60 | vfadd-vfdiv/vfsgnj*/vfmacc*/vmfeq*/vfsqrt/vfclass/vfcvt*/vfrsqrt7/vfrec7 |
| FP widening/narrowing | 33 | vfwadd/vfwsub/vfwmul/vfwmacc*/vfwcvt*/vfncvt* |
| Reductions | 16 | vredsum/vredmax/vfredosum/vfwredosum/... |
| Mask operations | 15 | vmand*/vcpop/vfirst/vid/viota/vmsbf/vmsif/vmsof |
| Permutation | 29 | vslide*/vrgather*/vcompress/vmv*/vfmv*/vmv*r |
| Edge cases | 24 | vl=0 no-op, tail undisturbed (full-VLMAX), LMUL>1 boundary, vle32ff fault, vill trap, fractional LMUL, tail/mask agnostic, stride zero/negative, ordered scatter, vxsat sticky, narrowing tail, widening LMUL>1, fflags |

## Building and running

### Prerequisites

- Python 3.10+ (for test generation)
- RISC-V GCC cross-compiler (`riscv64-unknown-linux-gnu-gcc`)
- A RISC-V machine with V extension (for running tests)

### Generate tests

```bash
python3 generate_tests.py -v
```

This also cleans the `tests/` directory before regenerating, so stale tests from previous runs are never left behind.

### Cross-compile

```bash
make build CC=riscv64-unknown-linux-gnu-gcc
```

### Run tests

```bash
bash scripts/run_tests.sh local
```

A nonzero exit code from any test binary identifies which check failed. Each `.S` file has comments at the top documenting what each check number means.

## Design

- **VLEN-agnostic**: Tests use `vsetvli` with small fixed vl (4 elements) and LMUL=1, working on any VLEN >= 128
- **No libc**: Each test is a standalone static binary using Linux syscall 93 (exit) directly
- **Generated**: A Python code generator produces all `.S` files to avoid duplication and ensure consistency
- **Side-effect checking**: After each instruction under test, CSRs and witness registers are verified
- **Full-VLMAX tail checking**: Edge case tests fill registers at VLMAX, dump with `vs1r.v`/`vs2r.v`, and verify tail bytes dynamically using `csrr vlenb`
- **vill/illegal-instruction testing**: Uses `clone`/`wait4` to fork a child process that executes a vector instruction with vill set, verifying the kernel kills it with SIGILL

## License

[CC0 1.0 Universal](LICENSE) - Public Domain
