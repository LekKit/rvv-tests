# Makefile for RVV 1.0 Test Suite
#
# Usage:
#   make                    - Generate + build all tests (native on riscv64)
#   make generate           - Generate tests from Python generator
#   make build              - Build all tests
#   make test               - Build and run all tests
#   make deploy             - Copy source tree to remote target
#   make remote-test        - Deploy, build, and run on remote target
#   make clean              - Remove build artifacts
#
# Configuration:
#   CC        - Compiler (default: gcc for native, riscv64-linux-gnu-gcc for cross)
#   ARCH      - Architecture string (default: rv64gcv)
#   ABI       - ABI (default: lp64d)
#   SSH_HOST  - SSH target for remote execution
#   SSH_PORT  - SSH port (default: 22)
#   REMOTE_DIR - Remote directory for test files

# ---- Defaults ----
ARCH       ?= rv64gcv
ABI        ?= lp64d
SSH_HOST   ?= kamillaova@lekkit.servebeer.com
SSH_PORT   ?= 2322
REMOTE_DIR ?= ~/rvv-tests
PYTHON     ?= python3

# ---- Detect if we're on RISC-V (native) or need cross-compilation ----
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_M),riscv64)
  CC      ?= gcc
  CFLAGS  = -march=$(ARCH) -mabi=$(ABI) -nostdlib -static -I include
else
  # Cross-compilation - try common prefixes
  CC      ?= riscv64-linux-gnu-gcc
  CFLAGS  = -march=$(ARCH) -mabi=$(ABI) -nostdlib -static -I include
endif

# ---- Source and build directories ----
INCLUDE_DIR = include
TEST_DIRS   = $(shell find tests -type d 2>/dev/null)
TEST_SRCS   = $(shell find tests -name '*.S' 2>/dev/null | sort)
TEST_BINS   = $(TEST_SRCS:.S=.bin)

# ---- Phony targets ----
.PHONY: all generate build test clean deploy remote-test remote-build \
        remote-run list help

all: generate build

help:
	@echo "RVV 1.0 Test Suite"
	@echo ""
	@echo "Targets:"
	@echo "  generate     - Generate test .S files from Python"
	@echo "  build        - Compile all .S files to binaries"
	@echo "  test         - Build and run all tests locally"
	@echo "  deploy       - rsync project to remote target"
	@echo "  remote-test  - Deploy + build + run on remote"
	@echo "  remote-build - Deploy + build on remote"
	@echo "  remote-run   - Run tests on remote (assumes built)"
	@echo "  list         - List all test files"
	@echo "  clean        - Remove build artifacts"
	@echo ""
	@echo "Configuration:"
	@echo "  SSH_HOST=$(SSH_HOST)"
	@echo "  SSH_PORT=$(SSH_PORT)"
	@echo "  CC=$(CC)"

# ---- Generate tests ----
generate:
	$(PYTHON) generate_tests.py -v

# ---- Build ----
# Each .S -> .bin
tests/%.bin: tests/%.S $(wildcard include/*.h)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ $<

build: $(TEST_BINS)

# ---- Run tests locally (must be on riscv64) ----
test: build
	@./scripts/run_tests.sh local

# ---- List tests ----
list:
	@echo "Test sources ($(words $(TEST_SRCS)) files):"
	@for f in $(TEST_SRCS); do echo "  $$f"; done

# ---- Remote operations ----
SSH_CMD = ssh -p $(SSH_PORT) $(SSH_HOST)
RSYNC_CMD = rsync -avz --delete -e "ssh -p $(SSH_PORT)"

deploy:
	$(RSYNC_CMD) \
		--exclude='*.bin' \
		--exclude='*.o' \
		--exclude='__pycache__' \
		--exclude='.git' \
		. $(SSH_HOST):$(REMOTE_DIR)/

remote-build: deploy
	$(SSH_CMD) "cd $(REMOTE_DIR) && make build"

remote-test: deploy
	$(SSH_CMD) "cd $(REMOTE_DIR) && make build && make test"

remote-run:
	$(SSH_CMD) "cd $(REMOTE_DIR) && ./scripts/run_tests.sh local"

# ---- Clean ----
clean:
	find tests -name '*.bin' -delete 2>/dev/null || true
	find tests -name '*.o' -delete 2>/dev/null || true
	rm -rf tests/generated 2>/dev/null || true
