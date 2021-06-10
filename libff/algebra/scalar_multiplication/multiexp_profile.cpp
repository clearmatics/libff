#include "libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp"
#include "libff/algebra/scalar_multiplication/multiexp.hpp"
#include "libff/common/profiling.hpp"
#include "libff/common/rng.hpp"

#include <cstdio>
#include <vector>

const size_t NUM_ITERATIONS = 10;
const size_t NUM_DIFFERENT_ELEMENTS = 32;

using namespace libff;

template<typename GroupT> using run_result_t = std::pair<long long, GroupT>;

template<typename T> using test_instances_t = std::vector<T>;

template<typename GroupT>
test_instances_t<GroupT> generate_group_elements(size_t num_elements)
{
    test_instances_t<GroupT> result;
    result.reserve(num_elements);
    assert(result.size() == 0);

    // Generating a random group element is expensive, so for now we only
    // generate NUM_DIFFERENT_ELEMENTS, and repeat them. Note, some methods
    // require input to be in special form.

    size_t i;
    for (i = 0; i < NUM_DIFFERENT_ELEMENTS; ++i) {
        GroupT x = GroupT::random_element();
        x.to_special();
        result.push_back(x);
    }
    assert(result.size() == NUM_DIFFERENT_ELEMENTS);

    for (; i < num_elements; ++i) {
        assert(result.size() == i);
        result.push_back(result[i % NUM_DIFFERENT_ELEMENTS]);
    }

    assert(result.size() == num_elements);
    return result;
}

template<typename FieldT>
test_instances_t<FieldT> generate_scalars(size_t num_elements)
{
    // Use SHA512_rng because it is much faster than FieldT::random_element()
    test_instances_t<FieldT> result;
    result.reserve(num_elements);
    for (size_t i = 0; i < num_elements; i++) {
        result.push_back(SHA512_rng<FieldT>(i));
    }

    assert(result.size() == num_elements);
    return result;
}

template<typename GroupT, typename FieldT, multi_exp_method Method>
run_result_t<GroupT> profile_multiexp(
    test_instances_t<GroupT> group_elements, test_instances_t<FieldT> scalars)
{
    long long start_time = get_nsec_time();

    GroupT answer;
    for (size_t iter = 0; iter < NUM_ITERATIONS; ++iter) {
        answer = multi_exp<GroupT, FieldT, Method>(
            group_elements.cbegin(),
            group_elements.cend(),
            scalars.cbegin(),
            scalars.cend(),
            1);
    }

    long long time_delta = get_nsec_time() - start_time;

    return run_result_t<GroupT>(time_delta, answer);
}

template<typename GroupT, typename FieldT>
void print_performance_csv(
    size_t expn_start,
    size_t expn_end_fast,
    size_t expn_end_naive,
    bool compare_answers)
{
    printf(
        "\t%16s\t%16s\t%16s\t%16s\t%16s\n",
        "bos-coster",
        "djb",
        "djb_signed",
        "djb_signed_mixed",
        "naive");
    for (size_t expn = expn_start; expn <= expn_end_fast; expn++) {
        printf("%ld", expn);
        fflush(stdout);

        test_instances_t<GroupT> group_elements =
            generate_group_elements<GroupT>(1 << expn);
        test_instances_t<FieldT> scalars = generate_scalars<FieldT>(1 << expn);

        run_result_t<GroupT> result_bos_coster =
            profile_multiexp<GroupT, FieldT, multi_exp_method_bos_coster>(
                group_elements, scalars);
        printf("\t%16lld", result_bos_coster.first);
        fflush(stdout);

        run_result_t<GroupT> result_djb =
            profile_multiexp<GroupT, FieldT, multi_exp_method_BDLO12>(
                group_elements, scalars);
        printf("\t%16lld", result_djb.first);
        fflush(stdout);

        if (compare_answers &&
            (result_bos_coster.second != result_djb.second)) {
            fprintf(stderr, "Answers NOT MATCHING (bos coster != djb)\n");
        }

        run_result_t<GroupT> result_djb_signed =
            profile_multiexp<GroupT, FieldT, multi_exp_method_BDLO12_signed>(
                group_elements, scalars);
        printf("\t%16lld", result_djb_signed.first);
        fflush(stdout);

        if (compare_answers &&
            (result_djb.second != result_djb_signed.second)) {
            fprintf(stderr, "Answers NOT MATCHING (djb != djb_signed)\n");
        }

        run_result_t<GroupT> result_djb_signed_mixed = profile_multiexp<
            GroupT,
            FieldT,
            multi_exp_method_BDLO12_signed_mixed>(group_elements, scalars);
        printf("\t%16lld", result_djb_signed_mixed.first);
        fflush(stdout);

        if (compare_answers &&
            (result_djb_signed.second != result_djb_signed_mixed.second)) {
            fprintf(
                stderr,
                "Answers NOT MATCHING (djb_signed != djb_signed_mixed)\n");
        }

        if (expn <= expn_end_naive) {
            run_result_t<GroupT> result_naive =
                profile_multiexp<GroupT, FieldT, multi_exp_method_naive>(
                    group_elements, scalars);
            printf("\t%16lld", result_naive.first);
            fflush(stdout);

            if (compare_answers &&
                (result_bos_coster.second != result_naive.second)) {
                fprintf(stderr, "Answers NOT MATCHING (bos coster != naive)\n");
            }
        }

        printf("\n");
    }
}

int main(void)
{
    print_compilation_info();

    alt_bn128_pp::init_public_params();
    printf("Profiling alt_bn128_G1\n");
    print_performance_csv<G1<alt_bn128_pp>, Fr<alt_bn128_pp> >(8, 20, 14, true);

    printf("Profiling alt_bn128_G2\n");
    print_performance_csv<G2<alt_bn128_pp>, Fr<alt_bn128_pp> >(8, 20, 14, true);

    return 0;
}
