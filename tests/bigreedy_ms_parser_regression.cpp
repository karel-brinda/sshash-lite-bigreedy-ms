#include "../include/dictionary.hpp"
#include "../include/builder/util.hpp"
#include "../include/bit_vector_iterator.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "../include/builder/parse_file.hpp"

namespace {

void assert_intervals(std::string const& sequence, uint64_t k,
                      std::initializer_list<sshash::embedding_interval> expected) {
    auto intervals = sshash::masked_superstring_intervals(sequence, k);
    if (intervals.size() != expected.size()) {
        throw std::runtime_error("masked-superstring scanner returned the wrong number of intervals");
    }
    auto expected_it = expected.begin();
    for (auto const& interval : intervals) {
        if (interval.begin != expected_it->begin || interval.end != expected_it->end) {
            throw std::runtime_error("masked-superstring scanner returned the wrong interval");
        }
        ++expected_it;
    }
}

void expect_invalid_tail_mask(std::string const& sequence, uint64_t k) {
    try {
        (void)sshash::masked_superstring_intervals(sequence, k);
    } catch (std::runtime_error const& error) {
        if (std::string(error.what()).find("final k-1 positions") != std::string::npos) return;
        throw;
    }
    throw std::runtime_error("uppercase mask position after the last valid start was accepted");
}

int run() {
    constexpr uint64_t k = 5;

    /* All of these have lowercase final k-1 positions. */
    assert_intervals("Acgtcg", k, {{0, 5}});      // first valid position only: [0, k)
    assert_intervals("aCgtcg", k, {{1, 6}});      // last valid position only: [n-k, n)
    assert_intervals("ACgtcg", k, {{0, 6}});      // run includes the last valid position
    assert_intervals("AcGtacga", k, {{0, 5}, {2, 7}});  // one lowercase position between runs
    expect_invalid_tail_mask("aCgtcG", k);  // position 5 is after the last valid start, 1

    auto root = std::filesystem::temp_directory_path() / "sshash-lite-bigreedy-ms-last-valid-start";
    std::filesystem::create_directories(root);
    try {
        sshash::build_configuration configuration;
        configuration.k = k;
        configuration.m = 3;
        configuration.tmp_dirname = root.string();

        std::string const masked_superstring = "aCgtcg";
        auto intervals = sshash::masked_superstring_intervals(masked_superstring, k);
        sshash::parse_data data(root.string());
        sshash::parse_bigreedy_intervals(masked_superstring, intervals, data, configuration);
        data.minimizers.merge();
        sshash::bit_vector_iterator strings_it(data.strings.strings, 2);

        if (data.num_kmers != 1 || data.minimizers.num_minimizers() != 1 ||
            data.strings.pieces != std::vector<uint64_t>{6} ||
            strings_it.read(2 * k) != sshash::util::string_to_uint64_no_reverse("CGTCG", k)) {
            throw std::runtime_error(
                "last valid masked k-mer start was not parsed into CGTCG with piece end 6");
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
    } catch (std::exception const& error) {
        std::cerr << "bigreedy-ms-parser-regression: " << error.what() << std::endl;
        return 1;
    }
}
