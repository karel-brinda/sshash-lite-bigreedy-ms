# End-to-end integration test for sshash-lite-bigreedy-ms:
#
#   unitigs (repetition-free SPSS)
#     --[kmercamel -S]-->      BIGREEDY masked superstring
#     --[sshash-lite build]--> membership index (+ --check)
#     --[sshash-lite query]--> query report
#
# Requirements: kmercamel on PATH (available on Bioconda), e.g.
#
#   conda install -c conda-forge -c bioconda kmercamel
#
# or
#
#   pixi global install -c conda-forge -c bioconda kmercamel
#
# (the bioconda channel must be given explicitly), and the
# sshash-lite-bigreedy-ms binary built in <repo>/build (see README).
#
# Usage:
#
#   make            # run the whole pipeline
#   make clean      # remove all generated files
#
# Tunables: K, M, T (query threshold), INPUT, QUERIES, BIN, OUT.
#
# Note: the input must be a repetition-free SPSS (unitigs or simplitigs;
# not matchtigs, which repeat k-mers). For an arbitrary k-mer set, compute
# simplitigs first, e.g. with prophasm (see README).

K ?= 31
M ?= 15
T ?= 0.8

ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../..)
BIN ?= $(ROOT)/build/sshash-lite-bigreedy-ms
INPUT ?= $(ROOT)/data/ecoli_unitigs_k31.fa.gz
QUERIES ?= $(ROOT)/data/queries/SRR9873306_1.10K.fastq.gz
OUT ?= output

.PHONY: all test check-tools clean

all: test

test: $(OUT)/bigreedy.sshash
	$(BIN) query $(OUT)/bigreedy.sshash $(QUERIES) -t $(T)
	@echo "=== INTEGRATION TEST PASSED ==="

$(OUT)/bigreedy.msfa: $(INPUT) | check-tools
	@mkdir -p $(OUT)
	kmercamel compute -k $(K) -S -o $@ $(INPUT)

$(OUT)/bigreedy.sshash: $(OUT)/bigreedy.msfa
	@mkdir -p $(OUT)/tmp
	$(BIN) build $< $(K) $(M) --canonical-parsing -d $(OUT)/tmp --check -o $@

check-tools:
	@command -v kmercamel >/dev/null 2>&1 || { \
		echo "error: kmercamel not found on PATH; install it e.g. with"; \
		echo "       conda install -c conda-forge -c bioconda kmercamel"; \
		echo "       or"; \
		echo "       pixi global install -c conda-forge -c bioconda kmercamel"; \
		exit 1; }
	@test -x "$(BIN)" || { \
		echo "error: $(BIN) not found;"; \
		echo "       build the project first (see the README)"; \
		exit 1; }

clean:
	rm -rf $(OUT)
