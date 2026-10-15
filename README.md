[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7277145.svg)](https://doi.org/10.5281/zenodo.7277145)

SSHash-Lite
===========

A membership-only version of [SSHash](https://github.com/jermp/sshash). It works with files with duplicate k-mers too.

For instructions on how to compile the code,
build the dictionaries, and data format, please refer to the instructions [here](https://github.com/jermp/sshash).

### Example usage

The default input format is a KmerCamel-style masked superstring: a FASTA file with one header
line and one unwrapped DNA sequence, where uppercase bases encode represented k-mer starts and
lowercase bases encode unrepresented starts. Normal KmerCamel `compute -o` output is directly
usable, without any conversion:

	kmercamel compute -k 31 -o bigreedy.msfa simplitigs.fa

	./sshash-lite-bigreedy-ms build bigreedy.msfa 31 15 --canonical-parsing -d tmp -o bigreedy.sshash --check

SSHash-Lite consumes the case mask during construction and stores the DNA sequence plus the
derived physical run ends, not the mask. KmerCamel's default reverse-complement behavior
corresponds to SSHash-Lite `--canonical-parsing`; KmerCamel `-u` corresponds to SSHash-Lite
regular parsing.

Query with

	./sshash-lite-bigreedy-ms query bigreedy.sshash ../data/queries/SRR9873306_1.10K.fastq.gz -t 0.8

where a query read is considered as positive if at least 80% of the k-mers are positive.

### SPSS input (one sequence per line)

To index a plain set of sequences instead (e.g. BCALM unitigs without duplicates, or
[matchtigs](https://github.com/algbio/matchtigs) with duplicates), pass `--spss`:

	./sshash-lite-bigreedy-ms build ../data/ecoli_matchtigs_k31.fa.gz 31 15 --spss -d tmp --check -o ecoli-matchtigs.sshash

**NOTE**: With `--spss`, input files are expected to have **one DNA sequence per line**. If a sequence spans multiple lines (e.g., multi-fasta), the lines should be concatenated before indexing.

To do so, the script `script/concat_lines.py` can be used:

	python3 script/concat_lines.py <multi_fasta>.fa.gz > output

### BIGREEDY superstring with explicit intervals

To index a BIGREEDY superstring with external interval annotations, you can provide
`--intervals embeddings.tsv` (this disables the masked-superstring default). The FASTA input must
contain exactly one unwrapped DNA sequence `B`. The tab-separated intervals file has one
zero-based, half-open `begin<TAB>end` interval per original SPSS string. Intervals must be
sorted by `begin`, have length at least `k`, and overlap by fewer than `k` bases.

	./sshash-lite-bigreedy-ms build B.fa 31 15 --intervals embeddings.tsv -d tmp -o B.sshash

Only k-mers fully contained in the listed intervals are indexed.

### Enable verbose output

By default we print a summary report for the whole query file.
To print the result for each query, define the flag `SSHASH_QUERY_VERBOSE_OUTPUT ` as follows:

	cmake .. -D SSHASH_QUERY_VERBOSE_OUTPUT=On
	make -j
    
and re-run the query.

The query output is written on `std::cerr` by default, so it is possible to capture the output with

	./sshash-lite-bigreedy-ms query ecoli-matchtigs.sshash ../data/queries/SRR9873306_1.10K.fastq.gz -t 0.8 2> query_output.txt
	
which writes the output to the file `query_output.txt`.
