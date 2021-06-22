#include "libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp"
#include "libff/algebra/curves/bls12_377/bls12_377_pp.hpp"
#include "libff/algebra/curves/curve_serialization.hpp"
#include "libff/common/profiling.hpp"

#include <aio.h>
#include <exception>
#include <fcntl.h>
#include <fstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>

using namespace libff;

static const size_t NUM_DIFFERENT_ELEMENTS = 1024;
static const size_t NUM_ELEMENTS_TO_READ = 1024 * 1024;
static const size_t NUM_ELEMENTS_IN_FILE = 256 * NUM_ELEMENTS_TO_READ;
static const size_t MAX_SPARSE_ELEMENT_INTERVAL =
    NUM_ELEMENTS_IN_FILE / NUM_ELEMENTS_TO_READ;

std::string get_filename(const std::string &identifier)
{
    return "group_elements_uncompressed_" + identifier + ".bin";
}

/// Returns true if the file was already present.
template<
    typename GroupT,
    form_t Form = form_montgomery,
    compression_t Comp = compression_off>
bool ensure_group_elements_file_uncompressed(const std::string &identifier)
{
    const std::string filename = get_filename(identifier);

    // If file doesn't exist, create it.
    struct stat s;
    if (stat(filename.c_str(), &s)) {
        std::cout << "  File '" << filename.c_str()
                  << "' does not exist. Creating ... ";
        std::flush(std::cout);

        // Fill a buffer with random elements
        std::vector<GroupT> elements;
        elements.reserve(NUM_DIFFERENT_ELEMENTS);
        for (size_t i = 0; i < NUM_DIFFERENT_ELEMENTS; ++i) {
            elements.push_back(GroupT::random_element());
        }

        // Use the buffer to fill the file
        std::ofstream out_s(
            filename.c_str(), std::ios_base::out | std::ios_base::binary);
        for (size_t i = 0; i < NUM_ELEMENTS_IN_FILE; ++i) {
            group_write<encoding_binary, Form, Comp>(
                elements[i % NUM_DIFFERENT_ELEMENTS], out_s);
        }
        out_s.close();

        std::cout << "Created\n";
        return false;
    }

    return true;
}

template<
    typename GroupT,
    form_t Form = form_montgomery,
    compression_t Comp = compression_off>
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
            for (size_t i = 0; i < NUM_ELEMENTS_TO_READ; ++i) {
                group_read<encoding_binary, Form, Comp>(
                    elements[i % NUM_DIFFERENT_ELEMENTS], in_s);
            }
            leave_block("Read group elements profiling");
        }

        in_s.close();
    }

    return true;
}

template<
    typename GroupT,
    form_t Form = form_montgomery,
    compression_t Comp = compression_off>
void profile_group_read_random_seek_uncompressed(const std::string &identifier)
{
    const std::string filename = get_filename(identifier);

    // Create a set of random indices. The i-th index to read from will be:
    //   (i + indices[i % NUM_DIFFERENT_ELEMENTS]) % NUM_ELEMENTS_IN_FILE
    std::vector<size_t> indices;
    indices.reserve(NUM_DIFFERENT_ELEMENTS);
    for (size_t i = 0; i < NUM_DIFFERENT_ELEMENTS; ++i) {
        indices.push_back(rand() % NUM_ELEMENTS_IN_FILE);
    }

    // Measure time taken to read the file
    std::cout << "  Random Access Read '" << filename.c_str() << "' ("
              << std::to_string(NUM_ELEMENTS_TO_READ) << " of "
              << std::to_string(NUM_ELEMENTS_IN_FILE) << " elements ...\n";
    {
        std::vector<GroupT> elements;
        elements.resize(NUM_DIFFERENT_ELEMENTS);

        std::ifstream in_s(
            filename.c_str(), std::ios_base::in | std::ios_base::binary);
        in_s.exceptions(
            std::ios_base::eofbit | std::ios_base::badbit |
            std::ios_base::failbit);

        const size_t element_size_on_disk = 2 * sizeof(elements[0].X);
        {
            enter_block("Read group elements profiling");
            for (size_t i = 0; i < NUM_ELEMENTS_TO_READ; ++i) {
                const size_t different_element_idx = i % NUM_DIFFERENT_ELEMENTS;
                const size_t idx =
                    (indices[different_element_idx] + i) % NUM_ELEMENTS_IN_FILE;
                const size_t offset = idx * element_size_on_disk;
                group_read<encoding_binary, Form, Comp>(
                    elements[different_element_idx], in_s.seekg(offset));
            }
            leave_block("Read group elements profiling");
        }

        in_s.close();
    }
}

template<
    typename GroupT,
    form_t Form = form_montgomery,
    compression_t Comp = compression_off>
void profile_group_read_random_seek_ordered_uncompressed(
    const std::string &identifier)
{
    const std::string filename = get_filename(identifier);

    // Measure time taken to read the file
    std::cout << "  Random Access Seeks (Ordered) '" << filename.c_str()
              << "' (" << std::to_string(NUM_ELEMENTS_TO_READ) << " of "
              << std::to_string(NUM_ELEMENTS_IN_FILE) << " elements ...\n";
    {
        std::vector<GroupT> elements;
        elements.resize(NUM_DIFFERENT_ELEMENTS);

        std::ifstream in_s(
            filename.c_str(), std::ios_base::in | std::ios_base::binary);
        in_s.exceptions(
            std::ios_base::eofbit | std::ios_base::badbit |
            std::ios_base::failbit);

        const size_t element_size_on_disk = 2 * sizeof(elements[0].X);
        const size_t element_interval_bytes =
            element_size_on_disk * NUM_ELEMENTS_IN_FILE / NUM_ELEMENTS_TO_READ;
        {
            enter_block("Read group elements profiling");
            for (size_t i = 0; i < NUM_ELEMENTS_TO_READ; ++i) {
                const size_t different_element_idx = i % NUM_DIFFERENT_ELEMENTS;
                const size_t offset = i * element_interval_bytes;
                group_read<encoding_binary, Form, Comp>(
                    elements[different_element_idx], in_s.seekg(offset));
            }
            leave_block("Read group elements profiling");
        }

        in_s.close();
    }
}

template<typename GroupT>
void profile_group_read_random_seek_fd_uncompressed(
    const std::string &identifier)
{
    const std::string filename = get_filename(identifier);

    // Create a set of random indices. The i-th index to read from will be:
    //   (i + indices[i % NUM_DIFFERENT_ELEMENTS]) % NUM_ELEMENTS_IN_FILE
    std::vector<size_t> indices;
    indices.reserve(NUM_DIFFERENT_ELEMENTS);
    for (size_t i = 0; i < NUM_DIFFERENT_ELEMENTS; ++i) {
        indices.push_back(rand() % NUM_ELEMENTS_IN_FILE);
    }

    // Measure time taken to read the file
    std::cout << "  Random Access Read '" << filename.c_str() << "' ("
              << std::to_string(NUM_ELEMENTS_TO_READ) << " of "
              << std::to_string(NUM_ELEMENTS_IN_FILE) << " elements ...\n";
    {
        std::vector<GroupT> elements;
        elements.resize(NUM_DIFFERENT_ELEMENTS);

        int f = open(filename.c_str(), O_RDONLY);
        if (f < 0) {
            throw std::runtime_error("failed to open " + filename);
        }

        const size_t element_size_on_disk = 2 * sizeof(elements[0].X);
        {
            enter_block("Read group elements profiling");
            for (size_t i = 0; i < NUM_ELEMENTS_TO_READ; ++i) {
                const size_t different_element_idx = i % NUM_DIFFERENT_ELEMENTS;
                const size_t idx =
                    (indices[different_element_idx] + i) % NUM_ELEMENTS_IN_FILE;
                const size_t offset = idx * element_size_on_disk;
                GroupT &dest_element = elements[different_element_idx];

                lseek(f, offset, SEEK_SET);
                read(f, &dest_element, element_size_on_disk);
            }
            leave_block("Read group elements profiling");
        }

        close(f);
    }
}

template<typename GroupT>
void profile_group_read_random_seek_fd_ordered_uncompressed(
    const std::string &identifier)
{
    const std::string filename = get_filename(identifier);

    // Measure time taken to read the file
    std::cout << "  Random Access Read '" << filename.c_str() << "' ("
              << std::to_string(NUM_ELEMENTS_TO_READ) << " of "
              << std::to_string(NUM_ELEMENTS_IN_FILE) << " elements ...\n";
    {
        std::vector<GroupT> elements;
        elements.resize(NUM_DIFFERENT_ELEMENTS);

        int f = open(filename.c_str(), O_RDONLY);
        if (f < 0) {
            throw std::runtime_error("failed to open " + filename);
        }

        const size_t element_size_on_disk = 2 * sizeof(elements[0].X);
        const size_t element_interval_bytes =
            element_size_on_disk * NUM_ELEMENTS_IN_FILE / NUM_ELEMENTS_TO_READ;
        {
            enter_block("Read group elements profiling");
            for (size_t i = 0; i < NUM_ELEMENTS_TO_READ; ++i) {
                const size_t different_element_idx = i % NUM_DIFFERENT_ELEMENTS;
                const size_t offset = i * element_interval_bytes;
                GroupT &dest_element = elements[different_element_idx];

                lseek(f, offset, SEEK_SET);
                read(f, &dest_element, element_size_on_disk);
            }
            leave_block("Read group elements profiling");
        }

        close(f);
    }
}

template<typename GroupT>
void profile_group_read_random_seek_mmap_ordered_uncompressed(
    const std::string &identifier)
{
    const std::string filename = get_filename(identifier);

    // Read sparse elements from a mem-mapped file, ensuring we always proceed
    // forwards.

    for (size_t interval = 1; interval < MAX_SPARSE_ELEMENT_INTERVAL;
         interval *= 2) {

        // Measure time taken to read the file
        std::cout << "  Random Access MMap (Ordered) '" << filename.c_str()
                  << "' (" << NUM_ELEMENTS_TO_READ << " elements of "
                  << (interval * NUM_ELEMENTS_TO_READ) << ") ...\n";
        {
            std::vector<GroupT> elements;
            elements.resize(NUM_DIFFERENT_ELEMENTS);

            const size_t element_size_on_disk = 2 * sizeof(elements[0].X);
            const size_t file_size =
                NUM_ELEMENTS_IN_FILE * element_size_on_disk;

            int fd = open(filename.c_str(), O_RDONLY);
            if (fd < 0) {
                throw std::runtime_error("failed to open " + filename);
            }

            const void *file_base = mmap(
                nullptr,
                file_size,
                PROT_READ,
                MAP_PRIVATE /* | MAP_NOCACHE */,
                fd,
                0);
            if (MAP_FAILED == file_base) {
                throw std::runtime_error(
                    std::string("mmap failed: ") + strerror(errno));
            }

            const size_t element_interval_bytes =
                element_size_on_disk * interval;

            {
                enter_block("Read group elements profiling");
                for (size_t i = 0; i < NUM_ELEMENTS_TO_READ; ++i) {
                    const size_t different_element_idx =
                        i % NUM_DIFFERENT_ELEMENTS;
                    const size_t offset = i * element_interval_bytes;
                    GroupT &dest_element = elements[different_element_idx];
                    const GroupT *src =
                        (const GroupT *)(const void
                                             *)((size_t)file_base + offset);

                    dest_element.X = src->X;
                    dest_element.Y = src->Y;
                }
                leave_block("Read group elements profiling");
            }

            if (0 != munmap((void *)file_base, file_size)) {
                throw std::runtime_error("munmap failed");
            }
            close(fd);
        }
    }
}

void cb_init(
    aiocb *cb, int fd, size_t offset_bytes, size_t size_bytes, void *dest)
{
    // std::cout << "cb_init: " << fd << ", off: " << offset_bytes << ", size: "
    // << size_bytes << ", dst: " << dest << "\n";

    memset(cb, 0, sizeof(aiocb));
    cb->aio_fildes = fd;
    cb->aio_offset = offset_bytes;
    cb->aio_buf = dest;
    cb->aio_nbytes = size_bytes;
    cb->aio_sigevent.sigev_notify = SIGEV_NONE;
    cb->aio_lio_opcode = LIO_READ;
}

void cb_enqueue(aiocb *cb)
{
    const int r = aio_read(cb);
    if (0 == r) {
        // std::cout << "cb_wait: enqueued\n";
        return;
    }

    throw std::runtime_error(
        std::string("error from aio_read: ") + strerror(errno));
}

ssize_t cb_wait(aiocb *cb)
{
    int err;
    for (;;) {
        err = aio_error(cb);
        if (EINPROGRESS == err) {
            std::this_thread::yield();
            // sleep(0);
            continue;
        }

        const ssize_t ret = aio_return(cb);
        if (0 <= ret) {
            // std::cout << "cb_wait: done\n";
            return ret;
        }

        if (ECANCELED == err) {
            throw std::runtime_error("aio_error: cancelled");
        }

        throw std::runtime_error(
            std::string("error from aio_error: ") + strerror(errno));
    }
}

template<typename GroupT>
void profile_group_read_random_aio_ordered_uncompressed(
    const std::string &identifier)
{
    const std::string filename = get_filename(identifier);

    // Treat the file as being divided into NUM_ELEMENTS_TO_READ groups, of
    // size NUM_ELEMENTS_IN_FILE / NUM_ELEMENTS_TO_READ, reading from a ranadom
    // offset within each group.

    const size_t GROUP_SIZE = NUM_ELEMENTS_IN_FILE / NUM_ELEMENTS_TO_READ;

    // Read from group_offset[i] in i-th group, i.e. global offset:
    //   GROUP_SIZE * i + group_offset[i]
    std::vector<size_t> group_offsets;
    group_offsets.reserve(NUM_ELEMENTS_TO_READ);
    for (size_t i = 0; i < NUM_ELEMENTS_TO_READ; ++i) {
        group_offsets.push_back((size_t)rand() * GROUP_SIZE / RAND_MAX);
    }

    int fd = open(filename.c_str(), O_RDONLY);
    if (fd < 0) {
        throw std::runtime_error("failed to open " + filename);
    }

    aiocb cb1;
    aiocb cb2;
    GroupT dest1;
    GroupT dest2;
    const size_t size_on_disk = 2 * sizeof(dest1.X);

    aiocb *cur_cb = &cb1;
    aiocb *next_cb = &cb2;
    GroupT *cur_dest = &dest1;
    GroupT *next_dest = &dest2;

    std::cout << "  Async Read '" << filename.c_str() << "' ("
              << std::to_string(NUM_ELEMENTS_TO_READ) << " of "
              << std::to_string(NUM_ELEMENTS_IN_FILE) << " elements ...\n";

    {
        enter_block("Read group elements profiling");

        // Enqueue first request
        cb_init(cur_cb, fd, group_offsets[0], size_on_disk, cur_dest);
        cb_enqueue(cur_cb);

        // Enqueue all requests
        for (size_t i = 1; i < NUM_ELEMENTS_TO_READ; ++i) {
            // Enqueue next (i-th)
            const size_t file_offset =
                (i * GROUP_SIZE + group_offsets[i]) * size_on_disk;
            cb_init(next_cb, fd, file_offset, size_on_disk, next_dest);
            cb_enqueue(next_cb);

            // Wait for current and process
            ssize_t bytes = cb_wait(cur_cb);
            if (bytes != size_on_disk) {
                throw std::runtime_error("insufficient bytes read");
            }
            // TODO: Check values

            // Swap (at which point, cur_cb is the next element to wait for, and
            // next_cb is unused).
            std::swap(cur_cb, next_cb);
            std::swap(cur_dest, next_dest);
        }

        // Wait for last request
        cb_wait(cur_cb);
        // TODO: Check values

        leave_block("Read group elements profiling");
    }
}

template<size_t BATCHSIZE> class batched_aio_reader
{
public:
    batched_aio_reader(int fd)
        : _fd(fd)
        , _active_batch(_batch1_ptrs)
        , _next_batch(_batch2_ptrs)
        , _next_batch_size(0)
    {
        memset(_batch1, 0, sizeof(_batch1));
        memset(_batch2, 0, sizeof(_batch2));
        for (size_t i = 0; i < BATCHSIZE; ++i) {
            _batch1_ptrs[i] = &_batch1[i];
            _batch2_ptrs[i] = &_batch2[i];
        }
    }

    void enqueue_read_first_batch(
        size_t offset_bytes, size_t size_bytes, void *dest)
    {
        assert(_next_batch_size < BATCHSIZE);

        aiocb *cb = _next_batch[_next_batch_size];
        cb_init(cb, _fd, offset_bytes, size_bytes, dest);

        // When first batch is full, enqueue it and allow writing to the next
        // batch.

        ++_next_batch_size;
        if (_next_batch_size == BATCHSIZE) {
            enqueue_next_batch();

            std::swap(_active_batch, _next_batch);
            _next_batch_size = 0;
        }
    }

    void enqueue_read(size_t offset_bytes, size_t size_bytes, void *dest)
    {
        assert(_next_batch_size < BATCHSIZE);

        aiocb *cb = _next_batch[_next_batch_size];
        cb_init(cb, _fd, offset_bytes, size_bytes, dest);

        ++_next_batch_size;
        if (_next_batch_size == BATCHSIZE) {
            enqueue_next_batch();
        }
    }

    void wait_last_read()
    {
        // Wait for the _active_batch
        for (size_t i = 0; i < BATCHSIZE; ++i) {
            int r = cb_wait(_active_batch[i]);
            if (0 > r) {
                throw std::runtime_error("bad read result");
            }
            if (0 == r) {
                std::cout << "+";
            }
        }

        // swap pointers so it can be written into
        std::swap(_active_batch, _next_batch);
        _next_batch_size = 0;
    }

protected:
    void enqueue_next_batch()
    {
        int r = lio_listio(LIO_NOWAIT, _next_batch, BATCHSIZE, nullptr);
        if (0 != r) {
            throw std::runtime_error("enqueue_batch error");
        }
    }

    const int _fd;
    aiocb _batch1[BATCHSIZE];
    aiocb _batch2[BATCHSIZE];
    aiocb *_batch1_ptrs[BATCHSIZE];
    aiocb *_batch2_ptrs[BATCHSIZE];
    aiocb **_active_batch; // Next batch to wait for
    aiocb **_next_batch;   // Batch being filled
    size_t _next_batch_size;
};

/// Similar to profile_group_read_random_aio_ordered_uncompressed, but with
/// more aio requests in-flight at the same time.
template<typename GroupT>
void profile_group_read_random_batch_aio_ordered_uncompressed(
    const std::string &identifier)
{
    const std::string filename = get_filename(identifier);

    // Treat the file as being divided into NUM_ELEMENTS_TO_READ groups, of
    // size NUM_ELEMENTS_IN_FILE / NUM_ELEMENTS_TO_READ, reading from a ranadom
    // offset within each group.

    const size_t GROUP_SIZE = NUM_ELEMENTS_IN_FILE / NUM_ELEMENTS_TO_READ;

    // Read from group_offset[i] in i-th group, i.e. global offset:
    //   GROUP_SIZE * i + group_offset[i]
    std::vector<size_t> group_offsets;
    group_offsets.reserve(NUM_ELEMENTS_TO_READ);
    for (size_t i = 0; i < NUM_ELEMENTS_TO_READ; ++i) {
        group_offsets.push_back((size_t)rand() * GROUP_SIZE / RAND_MAX);
    }

    int fd = open(filename.c_str(), O_RDONLY);
    if (fd < 0) {
        throw std::runtime_error("failed to open " + filename);
    }

    const size_t BATCH_SIZE = 8;
    static_assert(
        0 == (NUM_ELEMENTS_TO_READ % BATCH_SIZE), "invalid batch size");
    batched_aio_reader<BATCH_SIZE> reader(fd);

    GroupT dest1[BATCH_SIZE];
    GroupT dest2[BATCH_SIZE];
    const size_t size_on_disk = 2 * sizeof(dest1[0].X);

    GroupT *cur_dest = dest1;
    GroupT *next_dest = dest2;

    std::cout << "  Deep Async Read '" << filename.c_str() << "' ("
              << std::to_string(NUM_ELEMENTS_TO_READ) << " of "
              << std::to_string(NUM_ELEMENTS_IN_FILE) << " elements ...\n";

    {
        enter_block("Read group elements profiling");

        size_t i = 0;

        // Enqueue first requests
        for (size_t j = 0; j < BATCH_SIZE; ++j) {
            reader.enqueue_read_first_batch(
                group_offsets[i + j] * size_on_disk,
                size_on_disk,
                cur_dest + j);
        }
        i += BATCH_SIZE;

        // Enqueue all requests
        for (; i < NUM_ELEMENTS_TO_READ; i += BATCH_SIZE) {
            // Enqueue next batch
            for (size_t j = 0; j < BATCH_SIZE; ++j) {
                reader.enqueue_read(
                    ((i + j) * GROUP_SIZE + group_offsets[i + j]) *
                        size_on_disk,
                    size_on_disk,
                    cur_dest + j);
            }

            // Wait for current and process
            reader.wait_last_read();

            // Swap (at which point, cur_cb is the next element to wait for, and
            // next_cb is unused).
            std::swap(cur_dest, next_dest);
        }

        // Wait for last request
        reader.wait_last_read();

        leave_block("Read group elements profiling");
    }
}

template<typename GroupT> void run_profile(const std::string &identifier)
{
    std::cout << " profile: " << identifier << "\n";
    if (!ensure_group_elements_file_uncompressed<GroupT>(identifier)) {
        std::cout << "  Purge disk cache and re-run to profile.\n";
        return;
    }

    profile_group_read_sequential_uncompressed<GroupT>(identifier);
    // profile_group_read_random_seek_uncompressed<GroupT>(identifier);
    // profile_group_read_random_seek_ordered_uncompressed<GroupT>(identifier);
    // profile_group_read_random_seek_fd_uncompressed<GroupT>(identifier);
    // profile_group_read_random_seek_fd_ordered_uncompressed<GroupT>(identifier);
    profile_group_read_random_seek_mmap_ordered_uncompressed<GroupT>(
        identifier);
    // profile_group_read_random_aio_ordered_uncompressed<GroupT>(identifier);
    // profile_group_read_random_batch_aio_ordered_uncompressed<GroupT>(
    //     identifier);
}

int main(void)
{
    // Some configurations are disabled for now.

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
