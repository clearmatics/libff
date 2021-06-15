/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MULTIEXP_STREAM_TCC_
#define MULTIEXP_STREAM_TCC_

#include "libff/algebra/scalar_multiplication/multiexp.hpp"

namespace libff
{

template<form_t Form, compression_t Comp, typename GroupT>
void elements_from_stream_producer(
    std::istream &in_s,
    concurrent_fifo_spsc<GroupT> &fifo,
    const size_t num_entries)
{
    // Read and enqueue all entries from the stream
    for (size_t element_idx = 0; element_idx < num_entries; ++element_idx) {
        GroupT *dest = fifo.enqueue_begin_wait();
        group_read<encoding_binary, Form, Comp>(*dest, in_s);
        fifo.enqueue_end();
    }
}

template<typename GroupT, typename FieldT>
GroupT multi_exp_base_elements_from_fifo_all_rounds(
    concurrent_fifo_spsc<GroupT> &fifo,
    const std::vector<FieldT> &exponents,
    const size_t c)
{
    const size_t num_entries = exponents.size();
    const size_t num_digits = (FieldT::num_bits + c - 1) / c;
    const size_t num_buckets = 1 << (c - 1);

    // Allocate state for all rounds.
    std::vector<std::vector<GroupT>> round_buckets(num_digits);
    std::vector<std::vector<bool>> round_bucket_hit(num_digits);
    for (std::vector<GroupT> &buckets : round_buckets) {
        buckets.resize(num_buckets);
    }
    for (std::vector<bool> &bucket_hit : round_bucket_hit) {
        bucket_hit.resize(num_buckets);
        // bucket_hit.assign(num_buckets, false);
    }
    std::vector<ssize_t> digits(num_digits);

    // Process each element
    for (size_t el_idx = 0; el_idx < num_entries; ++el_idx) {
        // Decompose the scalar, and wait for an element from the fifo
        field_get_signed_digits(digits, exponents[el_idx], c, num_digits);
        const GroupT &group_element = *(fifo.dequeue_begin_wait());

        // Process all digits
        for (size_t digit_idx = 0; digit_idx < num_digits; ++digit_idx) {
            std::vector<GroupT> &buckets = round_buckets[digit_idx];
            std::vector<bool> &bucket_hit = round_bucket_hit[digit_idx];
            const ssize_t digit = digits[digit_idx];

            if (digit < 0) {
                const size_t bucket_idx = (-digit) - 1;
                assert(bucket_idx < num_buckets);
                if (bucket_hit[bucket_idx]) {
                    buckets[bucket_idx] =
                        buckets[bucket_idx].mixed_add(-group_element);
                } else {
                    buckets[bucket_idx] = -group_element;
                    bucket_hit[bucket_idx] = true;
                }
            } else if (digit > 0) {
                const size_t bucket_idx = digit - 1;
                assert(bucket_idx < num_buckets);
                if (bucket_hit[bucket_idx]) {
                    buckets[bucket_idx] =
                        buckets[bucket_idx].mixed_add(group_element);
                } else {
                    buckets[bucket_idx] = group_element;
                    bucket_hit[bucket_idx] = true;
                }
            }
        }

        fifo.dequeue_end();
    }

    // For each digit, sum the buckets and accumulate the total
    GroupT result = multiexp_accumulate_buckets(
        round_buckets[num_digits - 1],
        round_bucket_hit[num_digits - 1],
        num_buckets);
    round_buckets[num_digits - 1].clear();
    round_bucket_hit[num_digits - 1].clear();

    for (size_t digit_idx = num_digits - 2; digit_idx < num_digits;
         --digit_idx) {
        for (size_t i = 0; i < c; ++i) {
            result = result.dbl();
        }

        const GroupT digit_sum = multiexp_accumulate_buckets(
            round_buckets[digit_idx], round_bucket_hit[digit_idx], num_buckets);
        round_buckets[digit_idx].clear();
        round_bucket_hit[digit_idx].clear();

        result = result + digit_sum;
    }

    return result;
}

template<form_t Form, compression_t Comp, typename GroupT, typename FieldT>
GroupT multi_exp_stream(
    std::istream &base_elements_in, const std::vector<FieldT> &exponents)
{
    static const size_t FIFO_SIZE = 1024;
    const size_t num_entries = exponents.size();
    const size_t c = bdlo12_signed_optimal_c(num_entries);
    assert(c > 0);

    // Fifo for streaming
    concurrent_fifo_spsc<GroupT> fifo(FIFO_SIZE);

    // Launch the reading thread
    std::thread producer([&base_elements_in, &fifo, num_entries]() {
        elements_from_stream_producer<Form, Comp, GroupT>(
            base_elements_in, fifo, num_entries);
    });

    // Consume all elements from the fifo.
    const GroupT result =
        multi_exp_base_elements_from_fifo_all_rounds<GroupT, FieldT>(
            fifo, exponents, c);

    // Wait for reading thread
    producer.join();

    return result;
}

} // namespace libff

#endif // MULTIEXP_STREAM_TCC_
