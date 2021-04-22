#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/common/profiling.hpp>

#include <exception>
#include <fcntl.h>
#include <fstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace libff;

static const size_t NUM_DIFFERENT_ELEMENTS = 1024;
static const size_t NUM_ELEMENTS_TO_READ = 1024 * 1024;
static const size_t NUM_ELEMENTS_IN_FILE = 256 * NUM_ELEMENTS_TO_READ;

std::string get_filename(const std::string &identifier)
{
    return "group_elements_uncompressed_" + identifier + ".bin";
}

/// Returns true if the file was already present.
template<typename GroupT>
bool ensure_group_elements_file_uncompressed(const std::string &identifier)
{
    const std::string filename = get_filename(identifier);

    // If file doesn't exist, create it.
    struct stat s;
    if (stat(filename.c_str(), &s))
    {
        std::cout << "  File '" << filename.c_str()
                  << "' does not exist. Creating ... ";
        std::flush(std::cout);

        // Fill a buffer with random elements
        std::vector<GroupT> elements;
        elements.reserve(NUM_DIFFERENT_ELEMENTS);
        for (size_t i = 0; i < NUM_DIFFERENT_ELEMENTS; ++i)
        {
            elements.push_back(GroupT::random_element());
        }

        // Use the buffer to fill the file
        std::ofstream out_s(
            filename.c_str(), std::ios_base::out | std::ios_base::binary);
        for (size_t i = 0; i < NUM_ELEMENTS_IN_FILE; ++i)
        {
            elements[i % NUM_DIFFERENT_ELEMENTS].write_uncompressed(out_s);
        }
        out_s.close();

        std::cout << "Created\n";
        return false;
    }

    return true;
}

template<typename GroupT>
bool profile_group_read_sequential_uncompressed(const std::string &identifier)
{
    const std::string filename = get_filename(identifier);

    // Measure time taken to read the file
    std::cout << "  Sequential read '" << filename.c_str() << "' (expecting "
              << std::to_string(NUM_ELEMENTS_TO_READ) << " elements ...\n";
    {
        std::vector<GroupT> elements;
        elements.resize(NUM_DIFFERENT_ELEMENTS);

        std::ifstream in_s(
            filename.c_str(), std::ios_base::in | std::ios_base::binary);
        in_s.exceptions(
            std::ios_base::eofbit | std::ios_base::badbit |
            std::ios_base::failbit);

        {
            enter_block("Read group elements profiling");
            for (size_t i = 0; i < NUM_ELEMENTS_TO_READ; ++i)
            {
                GroupT::read_uncompressed(
                    in_s, elements[i % NUM_DIFFERENT_ELEMENTS]);
            }
            leave_block("Read group elements profiling");
        }

        in_s.close();
    }

    return true;
}

template<typename GroupT>
void run_profile(const std::string &identifier)
{
    std::cout << " profile: " << identifier << "\n";
    if (!ensure_group_elements_file_uncompressed<GroupT>(identifier))
    {
        std::cout << "  Purge disk cache and re-run to profile.\n";
        return;
    }

    profile_group_read_sequential_uncompressed<GroupT>(identifier);
}

int main(void)
{
    std::cout << "alt_bn128_pp\n";
    alt_bn128_pp::init_public_params();
    run_profile<alt_bn128_G1>("alt_bn128_G1");
    run_profile<alt_bn128_G2>("alt_bn128_G2");

    std::cout << "bls12_377_pp\n";
    bls12_377_pp::init_public_params();
    run_profile<bls12_377_G1>("bls12_377_G1");
    run_profile<bls12_377_G2>("bls12_377_G2");

    return 0;
}
