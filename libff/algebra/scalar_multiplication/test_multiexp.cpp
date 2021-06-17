#include "libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp"
#include "libff/algebra/curves/bls12_377/bls12_377_pp.hpp"
#include "libff/algebra/scalar_multiplication/multiexp.hpp"

#include <gtest/gtest.h>

using namespace libff;

namespace
{

template<typename GroupT, multi_exp_base_form BaseForm>
void test_multiexp_accumulate_buckets_group(const size_t num_buckets)
{
    using FieldT = typename GroupT::scalar_field;

    // Prepare values
    std::vector<GroupT> values;
    values.reserve(num_buckets);
    std::vector<bool> value_hit;
    value_hit.reserve(num_buckets);
    for (size_t i = 0; i < num_buckets; ++i) {
        const bool hit = (i & 1) || rand() % 2;
        value_hit.push_back(hit);
        values.push_back(GroupT::random_element());
    }

    // Expected value
    GroupT expected = GroupT::zero();
    for (size_t i = 0; i < num_buckets; ++i) {
        if (value_hit[i]) {
            expected = expected + (FieldT(i + 1) * values[i]);
        }
    }

    if (BaseForm == multi_exp_base_form_special) {
        batch_to_special(values);
    }

    // Actual value
    const GroupT actual =
        internal::multiexp_accumulate_buckets<GroupT, BaseForm>(
            values, value_hit, num_buckets);

    ASSERT_EQ(expected, actual);
}

template<typename GroupT, multi_exp_base_form BaseForm>
void test_multiexp_signed_digits_round_group(const size_t digit_idx)
{
    using FieldT = typename GroupT::scalar_field;
    using BigIntT =
        typename std::decay<decltype(FieldT::one().as_bigint())>::type;

    const size_t c = 4;
    const size_t num_buckets = 8; // 1 << (4 - 1);;

    // Prepare values
    std::vector<GroupT> bases{
        FieldT(1) * GroupT::one(),
        FieldT(2) * GroupT::one(),
        FieldT(3) * GroupT::one(),
        FieldT(3) * GroupT::one(),
        FieldT(3) * GroupT::one(),
        FieldT(4) * GroupT::one(),
    };
    std::vector<BigIntT> exponents{
        FieldT(4 << (c * digit_idx)).as_bigint(),
        FieldT(3 << (c * digit_idx)).as_bigint(),
        FieldT(2 << (c * digit_idx)).as_bigint(),
        FieldT(256 << (c * digit_idx)).as_bigint(), // ignored (digit = 0)
        FieldT(256 << (c * digit_idx)).as_bigint(), // ignored (digit = 0)
        FieldT(1 << (c * digit_idx)).as_bigint(),
    };

    // Prepare state (use 4 bits, to recover the original scalar values).
    const size_t num_entries = bases.size();
    std::vector<GroupT> buckets(num_buckets);
    std::vector<bool> bucket_hit(num_buckets);

    // Expected result is:
    //     4 * [1] + 3 * [2] + 2 * [3] + + 0 * [3] + 0 * [3] + 1 * [4]
    //   = [4] + [6] + [6] + 4
    //   = [20]
    const GroupT expected = FieldT(20) * GroupT::one();

    // Actual value
    const GroupT actual = internal::multi_exp_implementation<
        GroupT,
        FieldT,
        multi_exp_method_BDLO12_signed,
        BaseForm>::
        signed_digits_round(
            bases.begin(),
            bases.end(),
            exponents.begin(),
            buckets,
            bucket_hit,
            num_entries,
            num_buckets,
            c,
            digit_idx);

    ASSERT_EQ(expected, actual);
}

template<typename GroupT, multi_exp_method Method, multi_exp_base_form BaseForm>
void test_multiexp_inner()
{
    using FieldT = typename GroupT::scalar_field;

    const size_t c = 4;

    // Prepare values
    std::vector<GroupT> bases{
        FieldT(1) * GroupT::one(),
        FieldT(2) * GroupT::one(),
        FieldT(3) * GroupT::one(),
        FieldT(4) * GroupT::one(),
        FieldT(1) * GroupT::one(),
        FieldT(2) * GroupT::one(),
        FieldT(3) * GroupT::one(),
        FieldT(4) * GroupT::one(),
    };
    for (auto &b : bases) {
        b.to_affine_coordinates();
    }

    std::vector<FieldT> exponents{
        FieldT(4 << (c * 0)),
        FieldT(3 << (c * 1)),
        FieldT(2 << (c * 2)),
        FieldT(1 << (c * 3)),
        FieldT(4 << (c * 0)),
        FieldT(3 << (c * 1)),
        FieldT(2 << (c * 2)),
        FieldT(1 << (c * 3)),
    };

    // Expected result is:
    //     4 * [1] + (3 << c) * [2] + (2 << (2*c)) * [3] + (1 << (3*c)) * [4]
    //   = [4] + [6] + [6] + 4
    //   = [20]
    const long expected_unencoded =
        2 * ((4 << 0) + ((3 << c) * 2) + ((2 << (2 * c)) * 3) +
             ((1 << (3 * c)) * 4));
    const GroupT expected = FieldT(expected_unencoded) * GroupT::one();

    // Actual value
    const GroupT actual =
        internal::multi_exp_implementation<GroupT, FieldT, Method, BaseForm>::
            multi_exp_inner(
                bases.begin(), bases.end(), exponents.begin(), exponents.end());

    ASSERT_EQ(expected, actual);
}

template<typename GroupT> void test_multiexp_accumulate_buckets()
{
    test_multiexp_accumulate_buckets_group<GroupT, multi_exp_base_form_normal>(
        4);
    test_multiexp_accumulate_buckets_group<GroupT, multi_exp_base_form_normal>(
        16);
    test_multiexp_accumulate_buckets_group<GroupT, multi_exp_base_form_normal>(
        128);
    test_multiexp_accumulate_buckets_group<GroupT, multi_exp_base_form_special>(
        4);
    test_multiexp_accumulate_buckets_group<GroupT, multi_exp_base_form_special>(
        16);
    test_multiexp_accumulate_buckets_group<GroupT, multi_exp_base_form_special>(
        128);
}

template<typename GroupT> void test_multiexp_signed_digits_round()
{
    test_multiexp_signed_digits_round_group<GroupT, multi_exp_base_form_normal>(
        0);
    test_multiexp_signed_digits_round_group<
        GroupT,
        multi_exp_base_form_special>(0);
    test_multiexp_signed_digits_round_group<GroupT, multi_exp_base_form_normal>(
        1);
    test_multiexp_signed_digits_round_group<
        GroupT,
        multi_exp_base_form_special>(1);
    test_multiexp_signed_digits_round_group<GroupT, multi_exp_base_form_normal>(
        2);
    test_multiexp_signed_digits_round_group<
        GroupT,
        multi_exp_base_form_special>(2);
    test_multiexp_signed_digits_round_group<GroupT, multi_exp_base_form_normal>(
        3);
    test_multiexp_signed_digits_round_group<
        GroupT,
        multi_exp_base_form_special>(3);
}

// template<typename GroupT, multi_exp_method Method, multi_exp_base_form
// baseForm>
template<typename GroupT, multi_exp_method Method, multi_exp_base_form BaseForm>
void test_multi_exp_config(size_t num_elements)
{
    using Field = typename GroupT::scalar_field;

    std::vector<GroupT> base_elements;
    base_elements.reserve(num_elements);
    std::vector<Field> scalars;
    scalars.reserve(num_elements);

    // Create values for sum [n]*[1]_G + [n-1]*[2]_G + ... + [1]*[n]_G
    // where [x]_G is encoding of x in G. The result should then be
    GroupT g = GroupT::one();
    size_t expect_scalar = 0;
    for (size_t i = 0; i < num_elements; ++i) {
        base_elements.push_back(g);
        g = g + GroupT::one();

        scalars.push_back(num_elements - i);

        expect_scalar += (i + 1) * (num_elements - i);
    }
    const GroupT expect = Field(expect_scalar) * GroupT::one();

    if (BaseForm == multi_exp_base_form_special) {
        batch_to_special(base_elements);
    }

    // Call multi_exp
    const GroupT result1 = multi_exp<GroupT, Field, Method, BaseForm>(
        base_elements.begin(),
        base_elements.end(),
        scalars.begin(),
        scalars.end(),
        1);
    const GroupT result2 = multi_exp<GroupT, Field, Method, BaseForm>(
        base_elements.begin(),
        base_elements.end(),
        scalars.begin(),
        scalars.end(),
        2);
    const GroupT result4 = multi_exp<GroupT, Field, Method, BaseForm>(
        base_elements.begin(),
        base_elements.end(),
        scalars.begin(),
        scalars.end(),
        4);

    ASSERT_EQ(expect, result1);
    ASSERT_EQ(expect, result2);
    ASSERT_EQ(expect, result4);
}

template<typename GroupT, multi_exp_method Method>
void test_multi_exp_group_method()
{
    test_multi_exp_config<GroupT, Method, multi_exp_base_form_normal>(4);
    test_multi_exp_config<GroupT, Method, multi_exp_base_form_normal>(256);
    test_multi_exp_config<GroupT, Method, multi_exp_base_form_special>(4);
    test_multi_exp_config<GroupT, Method, multi_exp_base_form_special>(256);
}

template<typename GroupT> void test_multi_exp()
{
    test_multi_exp_group_method<GroupT, multi_exp_method_naive>();
    test_multi_exp_group_method<GroupT, multi_exp_method_naive_plain>();
    test_multi_exp_group_method<GroupT, multi_exp_method_bos_coster>();
    test_multi_exp_group_method<GroupT, multi_exp_method_BDLO12>();
    test_multi_exp_group_method<GroupT, multi_exp_method_BDLO12_signed>();
}

TEST(MultiExpTest, TestMultiExpAccumulateBucketsAltBN128)
{
    test_multiexp_accumulate_buckets<alt_bn128_G1>();
    test_multiexp_accumulate_buckets<alt_bn128_G2>();
}

TEST(MultiExpTest, TestMultiExpAccumulateBucketsBLS12_377)
{
    test_multiexp_accumulate_buckets<bls12_377_G1>();
    test_multiexp_accumulate_buckets<bls12_377_G2>();
}

TEST(MultiExpTest, TestMultiExpSignedDigitsRoundAltBN128)
{
    test_multiexp_signed_digits_round<alt_bn128_G1>();
    test_multiexp_signed_digits_round<alt_bn128_G2>();
}

TEST(MultiExpTest, TestMultiExpSignedDigitsRoundBLS12_377)
{
    test_multiexp_signed_digits_round<bls12_377_G1>();
    test_multiexp_signed_digits_round<bls12_377_G2>();
}

TEST(MultiExpTest, TestMultiExpAltBN128)
{
    test_multi_exp<alt_bn128_G1>();
    test_multi_exp<alt_bn128_G2>();
}

TEST(MultiExpTest, TestMultiExpBLS12_377)
{
    test_multi_exp<bls12_377_G1>();
    test_multi_exp<bls12_377_G2>();
}

} // namespace

int main(int argc, char **argv)
{
    libff::alt_bn128_pp::init_public_params();
    libff::bls12_377_pp::init_public_params();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
