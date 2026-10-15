#include "../include/dictionary.hpp"
#include "../include/query/streaming_query.hpp"

#include <array>
#include <cctype>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void write_file(std::filesystem::path const& filename, std::string const& contents) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("cannot create test input " + filename.string());
    out << contents;
}

std::string kmer_from_number(uint64_t value, uint64_t k) {
    static constexpr std::array<char, 4> alphabet = {'A', 'C', 'G', 'T'};
    std::string kmer(k, 'A');
    for (uint64_t i = 0; i != k; ++i) {
        kmer[k - i - 1] = alphabet[value & 3];
        value >>= 2;
    }
    return kmer;
}

std::string reverse_complement(std::string const& kmer) {
    std::string result(kmer.size(), 'A');
    for (uint64_t i = 0; i != kmer.size(); ++i) {
        switch (kmer[kmer.size() - i - 1]) {
            case 'A': result[i] = 'T'; break;
            case 'C': result[i] = 'G'; break;
            case 'G': result[i] = 'C'; break;
            case 'T': result[i] = 'A'; break;
            default: throw std::runtime_error("invalid DNA k-mer in test");
        }
    }
    return result;
}

void assert_differential(sshash::dictionary const& ordinary, sshash::dictionary const& intervals,
                         sshash::dictionary const& masked_superstring, bool canonical) {
    constexpr uint64_t k = 5;
    for (uint64_t i = 0; i != (1ULL << (2 * k)); ++i) {
        auto kmer = kmer_from_number(i, k);
        bool ordinary_answer = ordinary.is_member(kmer.c_str(), false);
        bool intervals_answer = intervals.is_member(kmer.c_str(), false);
        bool masked_superstring_answer = masked_superstring.is_member(kmer.c_str(), false);
        if (ordinary_answer != intervals_answer || ordinary_answer != masked_superstring_answer) {
            throw std::runtime_error(std::string(canonical ? "canonical" : "regular") +
                                     " ordinary, intervals, and masked-superstring membership differ for " +
                                     kmer);
        }
    }
}

std::string masked_sequence(std::string sequence, std::initializer_list<uint64_t> starts) {
    for (char& base : sequence) base = static_cast<char>(std::tolower(base));
    for (uint64_t start : starts) sequence[start] = static_cast<char>(std::toupper(sequence[start]));
    return sequence;
}

void expect_build_failure(std::filesystem::path const& input,
                          sshash::build_configuration const& configuration,
                          std::string const& expected_message) {
    try {
        sshash::dictionary dictionary;
        dictionary.build(input.string(), configuration);
    } catch (std::runtime_error const& error) {
        if (std::string(error.what()).find(expected_message) != std::string::npos) return;
        throw std::runtime_error("unexpected validation error: " + std::string(error.what()));
    }
    throw std::runtime_error("masked-superstring validation did not reject invalid input");
}

void assert_expected_membership(sshash::dictionary const& dictionary, bool canonical) {
    for (auto const& kmer : {"ACGTT", "CGTTA", "TAGTC", "AGTCG"}) {
        if (!dictionary.is_member(kmer, false)) {
            throw std::runtime_error("index missed expected k-mer " + std::string(kmer));
        }
        if (canonical and !dictionary.is_member(reverse_complement(kmer).c_str(), false)) {
            throw std::runtime_error("canonical index missed reverse complement of " +
                                     std::string(kmer));
        }
    }
    for (auto const& kmer : {"GTTAG", "TTAGT"}) {
        if (dictionary.is_member(kmer, false)) {
            throw std::runtime_error("BIGREEDY index included masked k-mer " + std::string(kmer));
        }
    }
}

int run() {
    constexpr uint64_t k = 5;
    auto root = std::filesystem::temp_directory_path() /
                ("sshash-lite-bigreedy-ms-" +
                 std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root);

    try {
        auto ordinary_fasta = root / "ordinary.fa";
        auto bigreedy_fasta = root / "bigreedy.fa";
        auto masked_fasta = root / "bigreedy.msfa";
        auto intervals = root / "embeddings.tsv";
        auto queries = root / "queries.fa";
        auto intervals_index = root / "bigreedy.sshash";
        auto masked_index = root / "bigreedy-masked.sshash";
        auto canonical_masked_index = root / "bigreedy-masked-canonical.sshash";
        write_file(ordinary_fasta, ">s1\nACGTTA\n>s2\nTAGTCG\n");
        write_file(bigreedy_fasta, ">B\nACGTTAGTCG\n");
        write_file(masked_fasta, ">B\n" + masked_sequence("ACGTTAGTCG", {0, 1, 4, 5}) + "\n");
        write_file(intervals, "0\t6\n4\t10\n");
        write_file(queries, ">q1\nACGTT\n>q2\nAGTCG\n");

        sshash::build_configuration configuration;
        configuration.k = k;
        configuration.m = 3;
        configuration.seed = 1;
        configuration.tmp_dirname = root.string();

        sshash::dictionary ordinary;
        ordinary.build(ordinary_fasta.string(), configuration);

        configuration.intervals_filename = intervals.string();
        {
            sshash::dictionary bigreedy;
            bigreedy.build(bigreedy_fasta.string(), configuration);
            essentials::save(bigreedy, intervals_index.c_str());
        }
        sshash::dictionary reloaded_bigreedy;
        essentials::load(reloaded_bigreedy, intervals_index.c_str());

        configuration.intervals_filename.clear();
        configuration.masked_superstring = true;
        {
            sshash::dictionary masked_superstring;
            masked_superstring.build(masked_fasta.string(), configuration);
            essentials::save(masked_superstring, masked_index.c_str());
        }
        sshash::dictionary reloaded_masked_superstring;
        essentials::load(reloaded_masked_superstring, masked_index.c_str());

        assert_differential(ordinary, reloaded_bigreedy, reloaded_masked_superstring, false);
        assert_expected_membership(ordinary, false);
        assert_expected_membership(reloaded_bigreedy, false);
        assert_expected_membership(reloaded_masked_superstring, false);

        auto report = reloaded_bigreedy.streaming_query_from_file(queries.string(), 1.0);
        if (report.num_kmers != 2 or report.num_positive_kmers != 2 or
            report.num_positive_reads != 2) {
            throw std::runtime_error("BIGREEDY streaming query did not find both valid k-mers");
        }

        configuration.canonical_parsing = true;
        configuration.intervals_filename.clear();
        configuration.masked_superstring = false;
        sshash::dictionary canonical_ordinary;
        canonical_ordinary.build(ordinary_fasta.string(), configuration);
        configuration.intervals_filename = intervals.string();
        sshash::dictionary canonical_bigreedy;
        canonical_bigreedy.build(bigreedy_fasta.string(), configuration);

        configuration.intervals_filename.clear();
        configuration.masked_superstring = true;
        {
            sshash::dictionary canonical_masked_superstring;
            canonical_masked_superstring.build(masked_fasta.string(), configuration);
            essentials::save(canonical_masked_superstring, canonical_masked_index.c_str());
        }
        sshash::dictionary canonical_reloaded_masked_superstring;
        essentials::load(canonical_reloaded_masked_superstring, canonical_masked_index.c_str());

        assert_differential(canonical_ordinary, canonical_bigreedy,
                            canonical_reloaded_masked_superstring, true);
        assert_expected_membership(canonical_ordinary, true);
        assert_expected_membership(canonical_bigreedy, true);
        assert_expected_membership(canonical_reloaded_masked_superstring, true);
        for (auto const& kmer : {"ACGTT", "CGTTA", "TAGTC", "AGTCG"}) {
            auto kmer_rc = reverse_complement(kmer);
            if (canonical_ordinary.is_member(kmer_rc.c_str(), false) !=
                canonical_bigreedy.is_member(kmer_rc.c_str(), false)) {
                throw std::runtime_error("canonical reverse-complement membership differs for " +
                                         std::string(kmer));
            }
        }

        auto invalid_tail = root / "invalid-tail.msfa";
        auto no_starts = root / "no-starts.msfa";
        write_file(invalid_tail, ">B\nACGTTAGTcG\n");
        write_file(no_starts, ">B\nacgttagtcg\n");
        configuration.canonical_parsing = false;
        expect_build_failure(invalid_tail, configuration, "final k-1 positions");
        expect_build_failure(no_starts, configuration, "no represented k-mer starts");

        auto isolated = root / "isolated.msfa";
        write_file(isolated, ">B\n" + masked_sequence("ACGTTAC", {0, 2}) + "\n");
        sshash::dictionary isolated_masked_superstring;
        isolated_masked_superstring.build(isolated.string(), configuration);
        if (!isolated_masked_superstring.is_member("ACGTT", false) ||
            !isolated_masked_superstring.is_member("GTTAC", false) ||
            isolated_masked_superstring.is_member("CGTTA", false)) {
            throw std::runtime_error("isolated uppercase mask positions were not interpreted correctly");
        }
    } catch (...) {
        std::filesystem::remove_all(root);
        throw;
    }
    std::filesystem::remove_all(root);
    return 0;
}

}  // namespace

int main() {
    try {
        return run();
    } catch (std::exception const& e) {
        std::cerr << "bigreedy-ms-regression: " << e.what() << std::endl;
        return 1;
    }
}
