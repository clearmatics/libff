// If NDEBUG is defined, assert turns into nothing. No checks are made and a
// lot of warnings are generated.
#ifdef NDEBUG
# undef NDEBUG
#endif

#include <iostream>
#include <sstream>
#include <string>

#include <vector>
#include <numeric>
#include <chrono>

#include <fstream>

#include <libff/common/profiling.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/algebra/curves/bw6_761/bw6_761_pp.hpp>

using namespace libff;

template<typename ppT>
std::string pairing_bench(size_t data_set_size, std::string curve_str)
{
    G1<ppT> G1_one = G1<ppT>::one();
    G2<ppT> G2_one = G2<ppT>::one();

    // Prepare test data
    printf("Running pairing benchmarks:\n");
    std::vector<G1<ppT>> Ps;
    std::vector<G2<ppT>> Qs;
    for (size_t i = 0; i < data_set_size; i++){
        Ps.push_back((Fr<ppT>::random_element()) * G1_one);
        Qs.push_back((Fr<ppT>::random_element()) * G2_one);
    }

    // Pairing computation
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < data_set_size; i++){
        ppT::reduced_pairing(Ps[i], Qs[i]);
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end-start;
    std::stringstream ss;
    ss << curve_str <<","<< data_set_size << "," << diff.count() << std::endl;
    return ss.str();
}

// TODO: Do a scalar mult benchmark

int main(void)
{
    // Init all params of the benchmarked curves
    alt_bn128_pp::init_public_params();
    bls12_377_pp::init_public_params();
    bw6_761_pp::init_public_params();
    mnt4_pp::init_public_params();
    mnt6_pp::init_public_params();

    std::ofstream file;
    std::string res;

    file.open("benchmark_data.txt");
    // CSV header
    file << "curve,pairing_nb,exec_time_sec\n";

    // Number of repetition of each run
    size_t repetition = 5;
    for (size_t i = 0; i < repetition; i++){
        // ALT_BN128 pairings
        res = pairing_bench<alt_bn128_pp>(1, "alt_bn128");
        file << res;
        res = pairing_bench<alt_bn128_pp>(10, "alt_bn128");
        file << res;
        res = pairing_bench<alt_bn128_pp>(100, "alt_bn128");
        file << res;
        //res = pairing_bench<alt_bn128_pp>(1000, "alt_bn128");
        //file << res;
        //res = pairing_bench<alt_bn128_pp>(10000, "alt_bn128");
        //file << res;
        //res = pairing_bench<alt_bn128_pp>(100000, "alt_bn128");
        //file << res;
        //res = pairing_bench<alt_bn128_pp>(1000000, "alt_bn128");
        //file << res;

        // BLS12-377 pairings
        res = pairing_bench<bls12_377_pp>(1, "bls12_377");
        file << res;
        res = pairing_bench<bls12_377_pp>(10, "bls12_377");
        file << res;
        res = pairing_bench<bls12_377_pp>(100, "bls12_377");
        file << res;
        //res = pairing_bench<bls12_377_pp>(1000, "bls12_377");
        //file << res;
        //res = pairing_bench<bls12_377_pp>(10000, "bls12_377");
        //file << res;

        // BW6-761 pairings
        res = pairing_bench<bw6_761_pp>(1, "bw6_761");
        file << res;
        res = pairing_bench<bw6_761_pp>(10, "bw6_761");
        file << res;
        res = pairing_bench<bw6_761_pp>(100, "bw6_761");
        file << res;

        // MNT4 pairings
        res = pairing_bench<mnt4_pp>(1, "mnt4_298");
        file << res;
        res = pairing_bench<mnt4_pp>(10, "mnt4_298");
        file << res;
        res = pairing_bench<mnt4_pp>(100, "mnt4_298");
        file << res;

        // MNT6 pairings
        res = pairing_bench<mnt6_pp>(1, "mnt6_298");
        file << res;
        res = pairing_bench<mnt6_pp>(10, "mnt6_298");
        file << res;
        res = pairing_bench<mnt6_pp>(100, "mnt6_298");
        file << res;
    }
    file.close();
}
