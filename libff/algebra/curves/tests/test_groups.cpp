/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

// If NDEBUG is defined, assert turns into nothing. No checks are made and a
// lot of warnings are generated.
#ifdef NDEBUG
# undef NDEBUG
#endif

#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
#include <libff/common/profiling.hpp>
#ifdef CURVE_BN128
#include <libff/algebra/curves/bn128/bn128_pp.hpp>
#endif
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/algebra/curves/bw6_761/bw6_761_pp.hpp>

#include <sstream>

using namespace libff;

template<typename GroupT>
void test_mixed_add()
{
    GroupT base, el, result;

    base = GroupT::zero();
    el = GroupT::zero();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::zero();
    el = GroupT::random_element();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::random_element();
    el = GroupT::zero();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::random_element();
    el = GroupT::random_element();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::random_element();
    el = base;
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base.dbl());
}

template<typename GroupT>
void test_group()
{
    bigint<1> rand1 = bigint<1>("76749407");
    bigint<1> rand2 = bigint<1>("44410867");
    bigint<1> randsum = bigint<1>("121160274");

    GroupT zero = GroupT::zero();
    assert(zero == zero);
    GroupT one = GroupT::one();
    assert(one == one);
    GroupT two = bigint<1>(2l) * GroupT::one();
    assert(two == two);
    GroupT five = bigint<1>(5l) * GroupT::one();

    GroupT three = bigint<1>(3l) * GroupT::one();
    GroupT four = bigint<1>(4l) * GroupT::one();

    assert(two+five == three+four);

    GroupT a = GroupT::random_element();
    GroupT b = GroupT::random_element();

    assert(one != zero);
    assert(a != zero);
    assert(a != one);

    assert(b != zero);
    assert(b != one);

    assert(a.dbl() == a + a);
    assert(b.dbl() == b + b);
    assert(one.add(two) == three);
    assert(two.add(one) == three);
    assert(a + b == b + a);
    assert(a - a == zero);
    assert(a - b == a + (-b));
    assert(a - b == (-b) + a);

    // handle special cases
    assert(zero + (-a) == -a);
    assert(zero - a == -a);
    assert(a - zero == a);
    assert(a + zero == a);
    assert(zero + a == a);

    assert((a + b).dbl() == (a + b) + (b + a));
    assert(bigint<1>("2") * (a + b) == (a + b) + (b + a));

    assert((rand1 * a) + (rand2 * a) == (randsum * a));

    assert(GroupT::order() * a == zero);
    assert(GroupT::order() * one == zero);
    assert((GroupT::order() * a) - a != zero);
    assert((GroupT::order() * one) - one != zero);

    test_mixed_add<GroupT>();
}

template<typename GroupT>
void test_mul_by_q()
{
    GroupT a = GroupT::random_element();
    assert((GroupT::base_field_char()*a) == a.mul_by_q());
}

template<typename GroupT>
void test_untwist_frobenius_twist()
{
    const GroupT a = GroupT::random_element() + GroupT::random_element();
    const GroupT uft = a.untwist_frobenius_twist();
    const GroupT uft_2 = uft.untwist_frobenius_twist();

    // \pi^2(P) - [t] \pi(P) + P = zero
    const typename GroupT::scalar_field trace_frob("9586122913090633730");

    const GroupT z = uft_2 - (trace_frob * uft) + (GroupT::base_field::mod * a);
    assert(z == GroupT::zero());
}

template<typename GroupT>
void test_mul_by_cofactor()
{
    const GroupT a = GroupT::random_element();
    const GroupT a_h = GroupT::h*a;
    assert(a_h == a.mul_by_cofactor());
}

template<typename GroupT>
void test_output()
{
    GroupT g = GroupT::zero();

    for (size_t i = 0; i < 1000; ++i)
    {
        std::stringstream ss;
        ss << g;
        g.write_compressed(ss);
        g.write_uncompressed(ss);

        GroupT gg;
        GroupT gg_comp;
        GroupT gg_uncomp;
        ss >> gg;
        GroupT::read_compressed(ss, gg_comp);
        GroupT::read_uncompressed(ss, gg_uncomp);

        assert(g == gg);
        assert(g == gg_comp);
        assert(g == gg_uncomp);
        /* use a random point in next iteration */
        g = GroupT::random_element();
    }
}

template<typename GroupT>
void test_group_membership_valid()
{
    for (size_t i = 0 ; i < 1000 ; ++i)
    {
        GroupT g = GroupT::random_element();
        assert(g.is_in_safe_subgroup());
    }
}

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

template<typename GroupT>
void test_group_membership_invalid_g1(const typename GroupT::base_field &x)
{
    const typename  GroupT::base_field x_squared = x * x;
    const typename  GroupT::base_field x_cubed = x_squared * x;
    const typename  GroupT::base_field y_squared =
        x_cubed + (GroupT::coeff_a * x_squared) + GroupT::coeff_b;
    const typename GroupT::base_field y = y_squared.sqrt();
    const GroupT g1_invalid(x, y, GroupT::base_field::one());

    assert(g1_invalid.is_well_formed());
    assert(!g1_invalid.is_in_safe_subgroup());
}

template<typename GroupT>
void test_group_membership_invalid_g2(const typename GroupT::twist_field &x)
{
    const typename GroupT::twist_field x_squared = x * x;
    const typename GroupT::twist_field x_cubed = x_squared * x;
    const typename GroupT::twist_field y_squared =
        x_cubed + (GroupT::coeff_a * x_squared) + GroupT::coeff_b;
    const typename GroupT::twist_field y = y_squared.sqrt();
    const GroupT g2_invalid(x, y, GroupT::twist_field::one());

    assert(g2_invalid.is_well_formed());
    assert(!g2_invalid.is_in_safe_subgroup());
}

template<typename ppT>
void test_check_membership()
{
    test_group_membership_valid<G1<ppT>>();
    test_group_membership_valid<G2<ppT>>();
}

template<>
void test_check_membership<alt_bn128_pp>()
{
    test_group_membership_valid<alt_bn128_G1>();
    test_group_membership_valid<alt_bn128_G2>();
    test_group_membership_invalid_g2<alt_bn128_G2>(alt_bn128_Fq2::one());
}

template<>
void test_check_membership<bls12_377_pp>()
{
    test_group_membership_valid<bls12_377_G1>();
    test_group_membership_valid<bls12_377_G2>();
    test_group_membership_invalid_g1<bls12_377_G1>(bls12_377_Fq(3));
    test_group_membership_invalid_g2<bls12_377_G2>(
        bls12_377_Fq(3) * bls12_377_Fq2::one());
    assert(profile_group_membership<bls12_377_G1>());
    assert(profile_group_membership<bls12_377_G2>());
}

template<>
void test_check_membership<bw6_761_pp>()
{
    test_group_membership_valid<bw6_761_G1>();
    test_group_membership_valid<bw6_761_G2>();
    test_group_membership_invalid_g1<bw6_761_G1>(bw6_761_Fq(6));
    test_group_membership_invalid_g2<bw6_761_G2>(bw6_761_Fq(0));
}

void
test_bls12_377()
{
    const bls12_377_G1 g1 = bls12_377_G1::random_element();

    // Ensure sigma endomorphism results in multiplication by expected lambda.
    const bls12_377_G1 sigma_g1 = g1.sigma();
    assert(
        (bls12_377_Fr("91893752504881257701523279626832445440") * g1) ==
        sigma_g1);
}

int main(void)
{
    std::cout << "edwards_pp\n";
    edwards_pp::init_public_params();
    test_group<G1<edwards_pp> >();
    test_output<G1<edwards_pp> >();
    test_group<G2<edwards_pp> >();
    test_output<G2<edwards_pp> >();
    test_mul_by_q<G2<edwards_pp> >();

    std::cout << "mnt4_pp\n";
    mnt4_pp::init_public_params();
    test_group<G1<mnt4_pp> >();
    test_output<G1<mnt4_pp> >();
    test_group<G2<mnt4_pp> >();
    test_output<G2<mnt4_pp> >();
    test_mul_by_q<G2<mnt4_pp> >();
    test_check_membership<mnt4_pp>();
    test_mul_by_cofactor<G1<mnt4_pp>>();
    test_mul_by_cofactor<G2<mnt4_pp>>();

    std::cout << "mnt6_pp\n";
    mnt6_pp::init_public_params();
    test_group<G1<mnt6_pp> >();
    test_output<G1<mnt6_pp> >();
    test_group<G2<mnt6_pp> >();
    test_output<G2<mnt6_pp> >();
    test_mul_by_q<G2<mnt6_pp> >();
    test_check_membership<mnt6_pp>();
    test_mul_by_cofactor<G1<mnt6_pp>>();
    test_mul_by_cofactor<G2<mnt6_pp>>();

    std::cout << "alt_bn128_pp\n";
    alt_bn128_pp::init_public_params();
    test_group<G1<alt_bn128_pp> >();
    test_output<G1<alt_bn128_pp> >();
    test_group<G2<alt_bn128_pp> >();
    test_output<G2<alt_bn128_pp> >();
    test_mul_by_q<G2<alt_bn128_pp> >();
    test_check_membership<alt_bn128_pp>();
    test_mul_by_cofactor<G1<alt_bn128_pp>>();
    test_mul_by_cofactor<G2<alt_bn128_pp>>();

    // Make sure that added curves pass the libff tests
    std::cout << "bls12_377_pp\n";
    bls12_377_pp::init_public_params();
    test_bls12_377();
    test_group<G1<bls12_377_pp> >();
    test_output<G1<bls12_377_pp> >();
    test_group<G2<bls12_377_pp> >();
    test_output<G2<bls12_377_pp> >();
    test_mul_by_q<G2<bls12_377_pp> >();
    test_check_membership<bls12_377_pp>();
    test_untwist_frobenius_twist<G2<bls12_377_pp>>();
    test_mul_by_cofactor<G1<bls12_377_pp>>();
    test_mul_by_cofactor<G2<bls12_377_pp>>();

    std::cout << "bw6_761_pp\n";
    bw6_761_pp::init_public_params();
    test_group<G1<bw6_761_pp> >();
    test_output<G1<bw6_761_pp> >();
    test_group<G2<bw6_761_pp> >();
    test_output<G2<bw6_761_pp> >();
    test_mul_by_q<G2<bw6_761_pp> >();
    test_check_membership<bw6_761_pp>();
    test_mul_by_cofactor<G1<bw6_761_pp>>();
    test_mul_by_cofactor<G2<bw6_761_pp>>();

// BN128 has fancy dependencies so it may be disabled
#ifdef CURVE_BN128
    std::cout << "bn128_pp\n";
    bn128_pp::init_public_params();
    test_group<G1<bn128_pp> >();
    test_output<G1<bn128_pp> >();
    test_group<G2<bn128_pp> >();
    test_output<G2<bn128_pp> >();
    test_check_membership<bn128_pp>();
    test_mul_by_cofactor<G1<bn128_pp>>();
    test_mul_by_cofactor<G2<bn128_pp>>();
#endif
}
