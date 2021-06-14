/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <exception>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/common/profiling.hpp>

using namespace libff;

/// Number of random elements to generate (since this process is expensive)
static const size_t ADD_NUM_DIFFERENT_ELEMENTS = 1024;

/// Number of element to operate on in a single iteration.
static const size_t ADD_NUM_ELEMENTS = 1024 * 1024;

template<typename GroupT> bool profile_group_add()
{
    std::vector<GroupT> elements;
    {
        elements.reserve(ADD_NUM_ELEMENTS);
        size_t i = 0;
        for (; i < ADD_NUM_DIFFERENT_ELEMENTS; ++i) {
            elements.push_back(GroupT::random_element());
        }
        for (; i < ADD_NUM_ELEMENTS; ++i) {
            elements.push_back(elements[i % ADD_NUM_DIFFERENT_ELEMENTS]);
        }
    }

    std::cout << "    num elements: " << std::to_string(ADD_NUM_ELEMENTS)
              << "\n";

    size_t num_elements = 0;
    GroupT accum = GroupT::zero();
    enter_block("group add operation profiling");
    for (const GroupT &el : elements) {
        accum = accum.add(el);
        num_elements++;
    }
    leave_block("group add operation profiling");

    if (num_elements != ADD_NUM_ELEMENTS) {
        throw std::runtime_error("invalid number of elements seen");
    }

    return true;
}

template<typename GroupT> bool profile_group_mixed_add()
{
    std::vector<GroupT> elements;
    {
        elements.reserve(ADD_NUM_ELEMENTS);
        size_t i = 0;
        for (; i < ADD_NUM_DIFFERENT_ELEMENTS; ++i) {
            GroupT e = GroupT::random_element();
            e.to_affine_coordinates();
            elements.push_back(e);
        }
        for (; i < ADD_NUM_ELEMENTS; ++i) {
            elements.push_back(elements[i % ADD_NUM_DIFFERENT_ELEMENTS]);
        }
    }

    std::cout << "    num elements: " << std::to_string(ADD_NUM_ELEMENTS)
              << "\n";

    size_t num_elements = 0;
    GroupT accum = GroupT::one();
    enter_block("group mixed add operation profiling");
    for (const GroupT &el : elements) {
        accum = accum.mixed_add(el);
        num_elements++;
    }
    leave_block("group mixed add operation profiling");

    if (num_elements != ADD_NUM_ELEMENTS) {
        throw std::runtime_error("invalid number of elements seen");
    }

    return true;
}

template<typename GroupT> bool profile_group_membership()
{
    static const size_t NUM_ELEMENTS = 1000;

    // Measure the time taken to check membership of 1000 elements. (Note all
    // elements are in fact members of the group - we are not testing
    // correctness here).

    std::vector<GroupT> elements;
    elements.reserve(NUM_ELEMENTS);
    for (size_t i = 0; i < NUM_ELEMENTS; ++i) {
        elements.push_back(GroupT::random_element());
    }

    enter_block("group membership profiling");
    for (const GroupT &el : elements) {
        if (!el.is_in_safe_subgroup()) {
            return false;
        }
    }
    leave_block("group membership profiling");

    return true;
}

int main(void)
{
    std::cout << "alt_bn128_pp\n";
    alt_bn128_pp::init_public_params();

    std::cout << "  profile_group_add<alt_bn128_G1>:\n";
    if (!profile_group_add<alt_bn128_G1>()) {
        throw std::runtime_error("failed");
    }

    std::cout << "  profile_group_mixed_add<alt_bn128_G1>:\n";
    if (!profile_group_mixed_add<alt_bn128_G1>()) {
        throw std::runtime_error("failed");
    }

    std::cout << "  profile_group_add<alt_bn128_G2>:\n";
    if (!profile_group_add<alt_bn128_G2>()) {
        throw std::runtime_error("failed");
    }

    std::cout << "  profile_group_mixed_add<alt_bn128_G2>:\n";
    if (!profile_group_mixed_add<alt_bn128_G2>()) {
        throw std::runtime_error("failed");
    }

    std::cout << "bls12_377_pp\n";
    bls12_377_pp::init_public_params();

    std::cout << "  profile_group_add<bls12_377_G1>:\n";
    if (!profile_group_add<bls12_377_G1>()) {
        throw std::runtime_error("failed");
    }

    std::cout << "  profile_group_mixed_add<bls12_377_G1>:\n";
    if (!profile_group_mixed_add<bls12_377_G1>()) {
        throw std::runtime_error("failed");
    }

    std::cout << "  profile_group_add<bls12_377_G2>:\n";
    if (!profile_group_add<bls12_377_G2>()) {
        throw std::runtime_error("failed");
    }

    std::cout << "  profile_group_membership<bls12_377_G1>:\n";
    if (!profile_group_membership<bls12_377_G1>()) {
        throw std::runtime_error("failed");
    }

    std::cout << "  profile_group_membership<bls12_377_G2>:\n";
    if (!profile_group_membership<bls12_377_G2>()) {
        throw std::runtime_error("failed");
    }

    return 0;
}
