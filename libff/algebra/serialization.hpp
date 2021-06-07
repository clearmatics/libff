/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef __LIBFF_ALGEBRA_SERIALIZATION_HPP__
#define __LIBFF_ALGEBRA_SERIALIZATION_HPP__

#include <stdint.h>
#include <string>

namespace libff
{

/// Encodings for (de)serialization.
enum encoding_t : uint8_t
{
    encoding_binary = 0,        // big endian
    encoding_json = 1,
};

/// Encodings for (de)serialization.
enum form_t : uint8_t
{
    form_plain = 0,
    form_montgomery = 1,
};

/// Enable / disable compression in (de)serialization.
enum compression_t : uint8_t
{
    compression_off = 0,
    compression_on = 1,
};

void hex_to_bytes_reversed(const std::string &hex, void *dest, size_t bytes);

std::string bytes_to_hex_reversed(const void *bytes, size_t num_bytes);

} // namespace

#endif // __LIBFF_ALGEBRA_SERIALIZATION_HPP__
