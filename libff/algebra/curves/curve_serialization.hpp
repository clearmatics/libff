/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef __LIBFF_ALGEBRA_CURVES_CURVE_SERIALIZATION_HPP__
#define __LIBFF_ALGEBRA_CURVES_CURVE_SERIALIZATION_HPP__

#include "libff/algebra/fields/field_serialization.hpp"

#include <iostream>

namespace libff
{

/// In-memory decompression (where we wish to separate the streaming of
/// compressed data, and decompression into a full group element.
template<form_t Form, typename GroupT>
void group_decompress(
    GroupT &v,
    const uint8_t buffer[field_binary_size<group_coord_type<GroupT>>()]);

template<encoding_t Enc, form_t Form, compression_t Comp, typename GroupT>
void group_read(GroupT &v, std::istream &in_s);

template<encoding_t Enc, form_t Form, compression_t Comp, typename GroupT>
void group_write(const GroupT &v, std::ostream &out_s);

} // namespace libff

#include "libff/algebra/curves/curve_serialization.tcc"

#endif // __LIBFF_ALGEBRA_CURVES_CURVE_SERIALIZATION_HPP__
