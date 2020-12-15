/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/common/profiling.hpp>

#include <exception>

using namespace libff;

template<typename GroupT>
bool profile_group_membership()
{
    static const size_t NUM_ELEMENTS = 1000;

    // Measure the time taken to check membership of 1000 elements. (Note all
    // elements are in fact members of the group - we are not testing
    // correctness here).

    std::vector<GroupT> elements;
    elements.reserve(NUM_ELEMENTS);
    for (size_t i = 0 ; i < NUM_ELEMENTS ; ++i)
    {
        elements.push_back(GroupT::random_element());
    }

    enter_block("group membership profiling");
    for (const GroupT &el : elements)
    {
        if (!el.is_in_safe_subgroup())
        {
            return false;
        }
    }
    leave_block("group membership profiling");

    return true;
}

int main(void)
{
    std::cout << "bls12_377_pp\n";
    bls12_377_pp::init_public_params();

    std::cout << "profile_group_membership<bls12_377_G1>:\n";
    if (!profile_group_membership<bls12_377_G1>()) {
        throw std::runtime_error("failed");
    }

    std::cout << "profile_group_membership<bls12_377_G2>:\n";
    if (!profile_group_membership<bls12_377_G2>()) {
        throw std::runtime_error("failed");
    }
}
