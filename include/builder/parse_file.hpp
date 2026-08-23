#pragma once

#include "../gz/zip_stream.hpp"

#include <limits>
#include <sstream>

namespace sshash {

struct parse_data {
    parse_data(std::string const& tmp_dirname) : num_kmers(0), minimizers(tmp_dirname) {}
    uint64_t num_kmers;
    minimizers_tuples minimizers;
    compact_string_pool strings;
};

struct embedding_interval {
    uint64_t begin;
    uint64_t end;
};

std::vector<embedding_interval> parse_embedding_intervals(std::string const& filename,
                                                           uint64_t superstring_size, uint64_t k) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening intervals file '" + filename + "'");

    std::vector<embedding_interval> intervals;
    std::string line;
    uint64_t line_number = 0;
    while (std::getline(is, line)) {
        ++line_number;
        auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos or line[first] == '#') continue;

        std::istringstream row(line);
        uint64_t begin, end;
        std::string extra;
        if (!(row >> begin >> end) or (row >> extra)) {
            throw std::runtime_error("malformed interval at " + filename + ":" +
                                     std::to_string(line_number));
        }
        if (begin > superstring_size or end > superstring_size or begin > end or
            begin > std::numeric_limits<uint64_t>::max() - k or begin + k > end) {
            throw std::runtime_error("invalid interval at " + filename + ":" +
                                     std::to_string(line_number));
        }
        if (!intervals.empty()) {
            auto const& previous = intervals.back();
            if (begin < previous.begin) {
                throw std::runtime_error("intervals must be sorted by begin coordinate");
            }
            if (begin < previous.end and previous.end - begin >= k) {
                throw std::runtime_error("consecutive intervals overlap by at least k");
            }
        }
        intervals.push_back({begin, end});
    }
    if (intervals.empty()) throw std::runtime_error("intervals file contains no intervals");
    return intervals;
}

std::string read_bigreedy_superstring(std::istream& is) {
    std::string header, sequence, extra;
    if (!std::getline(is, header) or header.empty() or header.front() != '>') {
        throw std::runtime_error("BIGREEDY input must start with one FASTA header");
    }
    if (!std::getline(is, sequence) or sequence.empty()) {
        throw std::runtime_error("BIGREEDY input must contain one unwrapped DNA sequence");
    }
    if (std::getline(is, extra)) {
        throw std::runtime_error("BIGREEDY input must contain exactly one unwrapped FASTA record");
    }
    if (!util::is_valid(sequence.data(), sequence.size())) {
        throw std::runtime_error("BIGREEDY superstring contains invalid DNA characters");
    }
    return sequence;
}

void parse_bigreedy(std::istream& is, parse_data& data, build_configuration const& build_config) {
    uint64_t k = build_config.k;
    uint64_t m = build_config.m;
    uint64_t seed = build_config.seed;
    uint64_t max_num_kmers_in_super_kmer = k - m + 1;

    if (max_num_kmers_in_super_kmer >= (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8))) {
        throw std::runtime_error(
            "max_num_kmers_in_super_kmer " + std::to_string(max_num_kmers_in_super_kmer) +
            " does not fit into " + std::to_string(sizeof(num_kmers_in_super_kmer_uint_type) * 8) +
            " bits");
    }

    std::string superstring = read_bigreedy_superstring(is);
    auto intervals = parse_embedding_intervals(build_config.intervals_filename, superstring.size(), k);
    uint64_t num_super_kmers = 0;

    for (auto const& interval : intervals) {
        uint64_t begin = 0;
        uint64_t end = 0;
        uint64_t interval_size = interval.end - interval.begin;
        uint64_t prev_minimizer = constants::invalid_uint64;

        auto append_super_kmer = [&]() {
            if (prev_minimizer == constants::invalid_uint64 or begin == end) return;

            uint64_t num_kmers_in_super_kmer = end - begin;
            for (uint64_t block_begin = 0; block_begin < num_kmers_in_super_kmer;
                 block_begin += max_num_kmers_in_super_kmer) {
                uint64_t num_kmers_in_block = std::min<uint64_t>(
                    max_num_kmers_in_super_kmer, num_kmers_in_super_kmer - block_begin);
                data.minimizers.emplace_back(prev_minimizer, interval.begin + begin + block_begin,
                                             num_kmers_in_block);
                ++num_super_kmers;
            }
        };

        while (end != interval_size - k + 1) {
            char const* kmer = superstring.data() + interval.begin + end;
            uint64_t uint64_kmer = util::string_to_uint64_no_reverse(kmer, k);
            uint64_t minimizer = util::compute_minimizer(uint64_kmer, k, m, seed);
            if (build_config.canonical_parsing) {
                uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, k);
                uint64_t minimizer_rc = util::compute_minimizer(uint64_kmer_rc, k, m, seed);
                minimizer = std::min<uint64_t>(minimizer, minimizer_rc);
            }

            if (prev_minimizer == constants::invalid_uint64) prev_minimizer = minimizer;
            if (minimizer != prev_minimizer) {
                append_super_kmer();
                begin = end;
                prev_minimizer = minimizer;
            }
            ++data.num_kmers;
            ++end;
        }
        append_super_kmer();
    }

    data.minimizers.finalize();
    std::vector<uint64_t> ends;
    ends.reserve(intervals.size());
    for (auto const& interval : intervals) ends.push_back(interval.end);
    data.strings.build_from_string(superstring, num_super_kmers, std::move(ends));

    std::cout << "read BIGREEDY superstring with " << superstring.size() << " bases and "
              << intervals.size() << " intervals, " << data.num_kmers << " kmers" << std::endl;
    std::cout << "num_kmers " << data.num_kmers << std::endl;
    std::cout << "num_super_kmers " << data.strings.num_super_kmers() << std::endl;
    std::cout << "num_pieces " << data.strings.pieces.size() << std::endl;
}

void parse_file(std::istream& is, parse_data& data, build_configuration const& build_config) {
    uint64_t k = build_config.k;
    uint64_t m = build_config.m;
    uint64_t seed = build_config.seed;
    uint64_t max_num_kmers_in_super_kmer = k - m + 1;
    uint64_t block_size = 2 * k - m;  // max_num_kmers_in_super_kmer + k - 1

    if (max_num_kmers_in_super_kmer >= (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8))) {
        throw std::runtime_error(
            "max_num_kmers_in_super_kmer " + std::to_string(max_num_kmers_in_super_kmer) +
            " does not fit into " + std::to_string(sizeof(num_kmers_in_super_kmer_uint_type) * 8) +
            " bits");
    }

    /* fit into the wanted number of bits */
    assert(max_num_kmers_in_super_kmer < (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8)));

    compact_string_pool::builder builder(k);

    std::string sequence;
    uint64_t prev_minimizer = constants::invalid_uint64;

    uint64_t begin = 0;  // begin of parsed super_kmer in sequence
    uint64_t end = 0;    // end of parsed super_kmer in sequence
    uint64_t num_sequences = 0;
    uint64_t num_bases = 0;
    bool glue = false;

    auto append_super_kmer = [&]() {
        if (sequence.empty() or prev_minimizer == constants::invalid_uint64 or begin == end) return;

        assert(end > begin);
        char const* super_kmer = sequence.data() + begin;
        uint64_t size = (end - begin) + k - 1;
        assert(util::is_valid(super_kmer, size));

        /* if num_kmers_in_super_kmer > k - m + 1, then split the super_kmer into blocks */
        uint64_t num_kmers_in_super_kmer = end - begin;
        uint64_t num_blocks = num_kmers_in_super_kmer / max_num_kmers_in_super_kmer +
                              (num_kmers_in_super_kmer % max_num_kmers_in_super_kmer != 0);
        assert(num_blocks > 0);
        for (uint64_t i = 0; i != num_blocks; ++i) {
            uint64_t n = block_size;
            if (i == num_blocks - 1) n = size;
            uint64_t num_kmers_in_block = n - k + 1;
            assert(num_kmers_in_block <= max_num_kmers_in_super_kmer);
            data.minimizers.emplace_back(prev_minimizer, builder.offset, num_kmers_in_block);
            builder.append(super_kmer + i * max_num_kmers_in_super_kmer, n, glue);
            if (glue) {
                assert(data.minimizers.back().offset > k - 1);
                data.minimizers.back().offset -= k - 1;
            }
            size -= max_num_kmers_in_super_kmer;
            glue = true;
        }
    };

    while (!is.eof()) {
        std::getline(is, sequence);  // header sequence
        std::getline(is, sequence);  // DNA sequence
        if (sequence.size() < k) continue;

        if (++num_sequences % 100000 == 0) {
            std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
                      << data.num_kmers << " kmers" << std::endl;
        }

        begin = 0;
        end = 0;
        glue = false;  // start a new piece
        prev_minimizer = constants::invalid_uint64;
        num_bases += sequence.size();

        while (end != sequence.size() - k + 1) {
            char const* kmer = sequence.data() + end;
            assert(util::is_valid(kmer, k));
            uint64_t uint64_kmer = util::string_to_uint64_no_reverse(kmer, k);
            uint64_t minimizer = util::compute_minimizer(uint64_kmer, k, m, seed);

            if (build_config.canonical_parsing) {
                uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, k);
                uint64_t minimizer_rc = util::compute_minimizer(uint64_kmer_rc, k, m, seed);
                minimizer = std::min<uint64_t>(minimizer, minimizer_rc);
            }

            if (prev_minimizer == constants::invalid_uint64) prev_minimizer = minimizer;
            if (minimizer != prev_minimizer) {
                append_super_kmer();
                begin = end;
                prev_minimizer = minimizer;
                glue = true;
            }

            ++data.num_kmers;
            ++end;
        }

        append_super_kmer();
    }

    data.minimizers.finalize();
    builder.finalize();
    builder.build(data.strings);

    std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
              << data.num_kmers << " kmers" << std::endl;
    std::cout << "num_kmers " << data.num_kmers << std::endl;
    std::cout << "num_super_kmers " << data.strings.num_super_kmers() << std::endl;
    std::cout << "num_pieces " << data.strings.pieces.size() << " (+"
              << (2.0 * data.strings.pieces.size() * (k - 1)) / data.num_kmers << " [bits/kmer])"
              << std::endl;
    assert(data.strings.pieces.size() == num_sequences);
}

parse_data parse_file(std::string const& filename, build_configuration const& build_config) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    std::cout << "reading file '" << filename << "'..." << std::endl;
    parse_data data(build_config.tmp_dirname);
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        if (build_config.intervals_filename.empty()) {
            parse_file(zis, data, build_config);
        } else {
            parse_bigreedy(zis, data, build_config);
        }
    } else {
        if (build_config.intervals_filename.empty()) {
            parse_file(is, data, build_config);
        } else {
            parse_bigreedy(is, data, build_config);
        }
    }
    is.close();
    return data;
}

}  // namespace sshash
