/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "libff/common/concurrent_fifo.hpp"

#include <gtest/gtest.h>
#include <thread>

using namespace libff;

namespace
{

static void test_concurrent_fifo(
    const size_t queue_size, const size_t total_num_values)
{
    concurrent_fifo_spsc<size_t> fifo(queue_size);
    size_t producer_blocks = 0;
    size_t consumer_blocks = 0;

    // thread to enqueue values in sequence
    auto fill_fn = [&fifo, &producer_blocks, total_num_values]() {
        for (size_t i = 0; i < total_num_values; ++i) {
            size_t *dest = fifo.try_enqueue_begin();
            if (!dest) {
                ++producer_blocks;
                do {
                    std::this_thread::yield();
                    dest = fifo.try_enqueue_begin();
                } while (!dest);
            }

            *dest = i;
            fifo.enqueue_end();
        }
    };

    std::thread filler(fill_fn);

    for (size_t i = 0; i < total_num_values; ++i) {
        const size_t *src = fifo.try_dequeue_begin();
        if (!src) {
            ++consumer_blocks;
            do {
                std::this_thread::yield();
                src = fifo.try_dequeue_begin();
            } while (!src);
        }

        ASSERT_EQ(i, *src);
        fifo.dequeue_end();
    }

    filler.join();

    // Uncomment to print contention stats
    // std::cout << "producer_blocks: " << producer_blocks
    //           << ", consumer_blocks: " << consumer_blocks << "\n";
}

TEST(CommonTests, ConcurrentFifoTest)
{
    // Many values, high contention
    test_concurrent_fifo(2, 1024 * 1024);

    // Larger queue, high throughpu
    test_concurrent_fifo(32, 8 * 1024 * 1024);
}

} // namespace
