/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef __LIBFF_ALGEBRA_FIELDS_FIELD_SERIALIZATION_HPP__
#define __LIBFF_ALGEBRA_FIELDS_FIELD_SERIALIZATION_HPP__

#include "libff/algebra/fields/bigint.hpp"
#include "libff/algebra/serialization.hpp"

namespace libff
{

template<typename BigIntT>
void bigint_from_hex(BigIntT &v, const std::string &hex);

template<typename BigIntT>
std::string bigint_to_hex(const BigIntT &v, bool prefix = false);

template<typename BigIntT>
void bigint_from_dec(BigIntT &v, const std::string &dec);

template<typename BigIntT> std::string bigint_to_dec(const BigIntT &v);

template<typename FieldT> constexpr size_t field_binary_size();

template<
    encoding_t Enc = encoding_binary,
    form_t Form = form_plain,
    typename FieldT>
void field_read(FieldT &v, std::istream &in_s);

template<
    encoding_t Enc = encoding_binary,
    form_t Form = form_plain,
    typename FieldT>
void field_write(const FieldT &v, std::ostream &out_s);

/// Read a field element with flags (2 bits) embedded in it.
template<form_t Form, typename FieldT>
void field_read_with_flags(FieldT &v, mp_limb_t &flags, std::istream &in_s);

/// Similar to the stream version, but extract from an in-memory buffer
template<form_t Form, typename FieldT>
void field_read_with_flags(
    FieldT &v,
    mp_limb_t &flags,
    const uint8_t buffer[field_binary_size<FieldT>()]);

template<form_t Form, typename FieldT>
void field_write_with_flags(
    const FieldT &v, mp_limb_t flags, std::ostream &out_s);

} // namespace libff

#include "libff/algebra/fields/field_serialization.tcc"

#endif // __LIBFF_ALGEBRA_FIELDS_FIELD_SERIALIZATION_HPP__
