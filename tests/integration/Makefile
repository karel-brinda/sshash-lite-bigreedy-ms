# End-to-end integration test for sshash-lite-bigreedy-ms:
#
#   SPSS (unitigs)
#     --[prophasm]-->          simplitigs
#     --[kmercamel -S]-->      BIGREEDY masked superstring
#     --[sshash-lite build]--> membership index (+ --check)
#     --[sshash-lite query]--> query report
#
# Requirements: prophasm and kmercamel on PATH (both are on Bioconda), e.g.
#
#   conda install -c conda-forge -c bioconda prophasm kmercamel
#
# or
#
#   pixi global install prophasm kmercamel
#
# and the sshash-lite-bigreedy-ms binary built in <repo>/build (see README).
#
# Usage:
#
#   make            # run the whole pipeline
#   make clean      # remove all generated files
#
# Tunables: K, M, T (query threshold), INPUT, QUERIES, BIN, OUT.

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

$(OUT)/simplitigs.fa: $(INPUT) | check-tools
	@mkdir -p $(OUT)
	prophasm -k $(K) -i $(INPUT) -o $@

$(OUT)/bigreedy.msfa: $(OUT)/simplitigs.fa
	kmercamel compute -k $(K) -S -o $@ $<

$(OUT)/bigreedy.sshash: $(OUT)/bigreedy.msfa
	@mkdir -p $(OUT)/tmp
	$(BIN) build $< $(K) $(M) --canonical-parsing -d $(OUT)/tmp --check -o $@

check-tools:
	@command -v prophasm >/dev/null 2>&1 || { \
		echo "error: prophasm not found on PATH;"; \
		echo "       install it e.g. with 'conda install -c conda-forge -c bioconda prophasm'"; \
		exit 1; }
	@command -v kmercamel >/dev/null 2>&1 || { \
		echo "error: kmercamel not found on PATH;"; \
		echo "       install it e.g. with 'conda install -c conda-forge -c bioconda kmercamel'"; \
		exit 1; }
	@test -x "$(BIN)" || { \
		echo "error: $(BIN) not found;"; \
		echo "       build the project first (see the README)"; \
		exit 1; }

clean:
	rm -rf $(OUT)
