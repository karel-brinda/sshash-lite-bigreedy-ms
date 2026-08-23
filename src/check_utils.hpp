#pragma once

#include "../include/gz/zip_stream.hpp"
#include "../include/builder/bigreedy_ms.hpp"

#include <sstream>

namespace sshash {

void random_kmer(char* kmer, uint64_t k) {
    for (uint64_t i = 0; i != k; ++i) kmer[i] = "ACGT"[rand() % 4];
}

bool check_correctness_membership(std::istream& is, dictionary const& dict) {
    uint64_t k = dict.k();
    uint64_t n = dict.size();

    std::string kmer_str(k, 0);
    std::cout << "checking correctness of membership..." << std::endl;

    std::string line;
    uint64_t num_kmers = 0;

    while (!is.eof()) {
        std::getline(is, line);  // header sequence
        std::getline(is, line);  // DNA sequence
        if (line.size() < k) continue;

        for (uint64_t i = 0; i + k <= line.size(); ++i) {
            assert(util::is_valid(line.data() + i, k));
            uint64_t uint64_kmer = util::string_to_uint64_no_reverse(line.data() + i, k);

            if (num_kmers != 0 and num_kmers % 5000000 == 0) {
                std::cout << "checked " << num_kmers << " kmers" << std::endl;
            }
            if ((num_kmers & 1) == 0) {
                /* transform 50% of the kmers into their reverse complements */
                uint64_kmer = util::compute_reverse_complement(uint64_kmer, k);
            }

            util::uint64_to_string_no_reverse(uint64_kmer, kmer_str.data(), k);
            bool answer = dict.is_member(kmer_str.c_str());

            if (!answer) {
                std::cout << "ERROR: got false but expected true for kmer '" << kmer_str << "'"
                          << std::endl;
            }
            assert(answer);

            ++num_kmers;
        }
    }
    std::cout << "checked " << num_kmers << " kmers" << std::endl;

    std::cout << "EVERYTHING OK!" << std::endl;

    std::cout << "checking correctness of negative membership with random kmers..." << std::endl;
    uint64_t num_queries = std::min<uint64_t>(1000000, n);
    for (uint64_t i = 0; i != num_queries; ++i) {
        random_kmer(kmer_str.data(), k);
        /*
            We could use a std::unordered_set to check if kmer is really absent,
            but that would take much more memory...
        */
        bool answer = dict.is_member(kmer_str.c_str());
        if (answer) std::cerr << "kmer '" << kmer_str << "' found!" << std::endl;
    }

    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

/*
   The input file must be the one the index was built from.
*/
bool check_correctness_membership(dictionary const& dict, std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        good = check_correctness_membership(zis, dict);
    } else {
        good = check_correctness_membership(is, dict);
    }
    is.close();
    return good;
}

bool check_correctness_membership_bigreedy(std::istream& is, dictionary const& dict,
                                           std::string const& intervals_filename) {
    uint64_t k = dict.k();
    std::string header, superstring;
    if (!std::getline(is, header) or !std::getline(is, superstring)) {
        throw std::runtime_error("cannot read BIGREEDY superstring for checking");
    }

    std::ifstream intervals_file(intervals_filename.c_str());
    if (!intervals_file.good()) {
        throw std::runtime_error("cannot open intervals file '" + intervals_filename + "'");
    }

    std::cout << "checking BIGREEDY membership correctness..." << std::endl;
    std::string line;
    std::string kmer_str(k, 0);
    uint64_t num_kmers = 0;
    while (std::getline(intervals_file, line)) {
        auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos or line[first] == '#') continue;
        std::istringstream row(line);
        uint64_t begin, end;
        if (!(row >> begin >> end)) {
            throw std::runtime_error("malformed interval while checking BIGREEDY membership");
        }
        for (uint64_t pos = begin; pos + k <= end; ++pos) {
            char const* kmer = superstring.data() + pos;
            if (dict.canonicalized() and (num_kmers & 1) == 0) {
                uint64_t uint64_kmer = util::string_to_uint64_no_reverse(kmer, k);
                uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, k);
                util::uint64_to_string_no_reverse(uint64_kmer_rc, kmer_str.data(), k);
                kmer = kmer_str.data();
            }
            bool answer = dict.is_member(kmer, false);
            if (!answer) {
                throw std::runtime_error("BIGREEDY index missed a k-mer during correctness check");
            }
            ++num_kmers;
        }
    }
    std::cout << "checked " << num_kmers << " kmers" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

bool check_correctness_membership_masked_superstring(std::istream& is, dictionary const& dict) {
    uint64_t k = dict.k();
    std::string masked_superstring = read_masked_superstring(is);
    auto intervals = masked_superstring_intervals(masked_superstring, k);

    std::cout << "checking masked-superstring membership correctness..." << std::endl;
    std::string kmer_str(k, 0);
    uint64_t num_kmers = 0;
    for (auto const& interval : intervals) {
        for (uint64_t pos = interval.begin; pos + k <= interval.end; ++pos) {
            char const* kmer = masked_superstring.data() + pos;
            if (dict.canonicalized() and (num_kmers & 1) == 0) {
                uint64_t uint64_kmer = util::string_to_uint64_no_reverse(kmer, k);
                uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, k);
                util::uint64_to_string_no_reverse(uint64_kmer_rc, kmer_str.data(), k);
                kmer = kmer_str.data();
            }
            if (!dict.is_member(kmer, false)) {
                throw std::runtime_error(
                    "masked-superstring index missed a k-mer during correctness check");
            }
            ++num_kmers;
        }
    }
    std::cout << "checked " << num_kmers << " kmers" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

bool check_correctness_membership_masked_superstring(dictionary const& dict,
                                                      std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        good = check_correctness_membership_masked_superstring(zis, dict);
    } else {
        good = check_correctness_membership_masked_superstring(is, dict);
    }
    is.close();
    return good;
}

bool check_correctness_membership_bigreedy(dictionary const& dict, std::string const& filename,
                                           std::string const& intervals_filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        good = check_correctness_membership_bigreedy(zis, dict, intervals_filename);
    } else {
        good = check_correctness_membership_bigreedy(is, dict, intervals_filename);
    }
    is.close();
    return good;
}

}  // namespace sshash
