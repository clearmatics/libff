#include <gtest/gtest.h>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>

using namespace libff;

namespace
{

template<typename GroupT>
void test_multiexp_accumulate_buckets(const size_t num_buckets)
{
    using FieldT = typename GroupT::scalar_field;

    // Prepare values
    std::vector<GroupT> values;
    values.reserve(num_buckets);
    std::vector<bool> value_hit;
    value_hit.reserve(num_buckets);
    for (size_t i = 0; i < num_buckets; ++i) {
        const bool hit = rand() % 2;
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

    // Actual value
    const GroupT actual =
        multiexp_accumulate_buckets(values, value_hit, num_buckets);

    ASSERT_EQ(expected, actual);
}

template<typename GroupT, bool MixedAddition>
void test_multiexp_signed_digits_round(const size_t digit_idx)
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
    const GroupT actual =
        multiexp_signed_digits_round<GroupT, BigIntT, MixedAddition>(
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

template<typename GroupT, multi_exp_method Method> void test_multiexp_inner()
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
    const GroupT actual = multi_exp_inner<GroupT, FieldT, Method>(
        bases.begin(), bases.end(), exponents.begin(), exponents.end());

    ASSERT_EQ(expected, actual);
}

TEST(MultiExpTest, TestMultiExpAccumulateBuckets)
{
    test_multiexp_accumulate_buckets<alt_bn128_G1>(4);
    test_multiexp_accumulate_buckets<alt_bn128_G1>(16);
    test_multiexp_accumulate_buckets<alt_bn128_G1>(128);
}

TEST(MultiExpTest, TestMultiExpSignedDigitsRound)
{
    test_multiexp_signed_digits_round<alt_bn128_G1, false>(0);
    test_multiexp_signed_digits_round<alt_bn128_G1, true>(0);
    test_multiexp_signed_digits_round<alt_bn128_G1, false>(1);
    test_multiexp_signed_digits_round<alt_bn128_G1, true>(1);
    test_multiexp_signed_digits_round<alt_bn128_G1, false>(2);
    test_multiexp_signed_digits_round<alt_bn128_G1, true>(2);
    test_multiexp_signed_digits_round<alt_bn128_G1, false>(3);
    test_multiexp_signed_digits_round<alt_bn128_G1, true>(3);
}

TEST(MultiExpTest, TestMultiExpInner)
{
    test_multiexp_inner<alt_bn128_G1, multi_exp_method_BDLO12_signed_mixed>();
}

} // namespace

int main(int argc, char **argv)
{
    libff::alt_bn128_pp::init_public_params();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
