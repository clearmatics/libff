/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef __LIBFF_ALGEBRA_FIELDS_FIELD_SERIALIZATION_TCC__
#define __LIBFF_ALGEBRA_FIELDS_FIELD_SERIALIZATION_TCC__

#include "libff/algebra/serialization.hpp"

#include <algorithm>

namespace libff
{

namespace internal
{

template<encoding_t Enc, typename FieldT, form_t Form = form_plain>
class field_element_codec;

template<typename FieldT, form_t Form>
class field_element_codec<encoding_json, FieldT, Form>
{
public:
    // Convert a field element to JSON
    static void write(const FieldT &field_el, std::ostream &out_s)
    {
        using base_field_t =
            typename std::decay<decltype(field_el.coeffs[0])>::type;

        // Note that we write components of extension fields
        // highest-order-first.
        out_s << '[';
        size_t i = FieldT::tower_extension_degree - 1;
        do {
            // out_s << field_element_to_json(field_el.coeffs[i]);
            field_element_codec<encoding_json, base_field_t, Form>::write(
                field_el.coeffs[i], out_s);
            if (i > 0) {
                out_s << ',';
            }
        } while (i-- > 0);
        out_s << ']';
    }

    // Read a field element from JSON
    static void read(FieldT &field_el, std::istream &in_s)
    {
        using base_field_t =
            typename std::decay<decltype(field_el.coeffs[0])>::type;

        // Read opening '[' char, then each component (highest-order-first)
        // separated by ',' char, then a closing ']' char.

        char sep;
        in_s >> sep;
        if (sep != '[') {
            throw std::runtime_error("expected opening bracket");
        }

        size_t i = FieldT::tower_extension_degree - 1;
        do {
            field_element_codec<encoding_json, base_field_t, Form>::read(
                field_el.coeffs[i], in_s);
            // field_element_read_json(field_el.coeffs[i], in_s);
            if (i > 0) {
                in_s >> sep;
                if (sep != ',') {
                    throw std::runtime_error("expected comma separator");
                }
            }
        } while (i-- > 0);

        in_s >> sep;
        if (sep != ']') {
            throw std::runtime_error("expected closing bracket");
        }
    }
};

// Implementation of field_element_codec<encoding_json, ...> for the base-case
// of Fp_model types.
template<mp_size_t n, const libff::bigint<n> &modulus, form_t Form>
class field_element_codec<encoding_json, libff::Fp_model<n, modulus>, Form>
{
public:
    using Field = libff::Fp_model<n, modulus>;
    static void write(const Field &field_el, std::ostream &out_s)
    {
        if (Form == form_plain) {
            out_s << '"' << bigint_to_hex(field_el.as_bigint()) << '"';
        } else {
            out_s << '"' << bigint_to_hex(field_el.mont_repr) << '"';
        }
    };
    static void read(Field &field_el, std::istream &in_s)
    {
        char quote;
        in_s >> quote;
        if (quote != '"') {
            throw std::runtime_error("expected json string");
        }
        std::string bigint_hex;
        try {
            std::getline(in_s, bigint_hex, '"');
        } catch (...) {
            throw std::runtime_error("json string not terminated");
        }

        if (Form == form_plain) {
            bigint<n> plain;
            bigint_from_hex(plain, bigint_hex);
            field_el = Field(plain);
        } else {
            bigint_from_hex(field_el.mont_repr, bigint_hex);
        }
    }
};

// Generic reader and write for fields and field extensions.
template<typename FieldT, form_t Form>
class field_element_codec<encoding_binary, FieldT, Form>
{
public:
    static void write(const FieldT &field_el, std::ostream &out_s)
    {
        using base_field_t =
            typename std::decay<decltype(field_el.coeffs[0])>::type;
        for (size_t i = 0; i < FieldT::tower_extension_degree; ++i) {
            // field_element_write_bytes(field_el.coeffs[i], out_s);
            field_element_codec<encoding_binary, base_field_t, Form>::write(
                field_el.coeffs[i], out_s);
        }
    }
    static void read(FieldT &field_el, std::istream &in_s)
    {
        using base_field_t =
            typename std::decay<decltype(field_el.coeffs[0])>::type;
        for (size_t i = 0; i < FieldT::tower_extension_degree; ++i) {
            // field_element_read_bytes(field_el.coeffs[i], in_s);
            field_element_codec<encoding_binary, base_field_t, Form>::read(
                field_el.coeffs[i], in_s);
        }
    }
};

/// Implementation of field_element_bytes for the base-case of Fp_model types.
/// Big-endian bigint values (i.e. not in montgomery form).
template<mp_size_t n, const libff::bigint<n> &modulus, form_t Form>
class field_element_codec<encoding_binary, libff::Fp_model<n, modulus>, Form>
{
public:
    using Field = libff::Fp_model<n, modulus>;
    static void write(const Field &field_el, std::ostream &out_s)
    {
        // Convert to bigint, reverse bytes in-place, and write to stream.
        libff::bigint<n> bi;
        if (Form == form_plain) {
            bi = field_el.as_bigint();
        } else {
            bi = field_el.mont_repr;
        }
        std::reverse((char *)(&bi), (char *)(&bi + 1));
        out_s.write((const char *)(&bi.data[0]), sizeof(bi));
    }
    static void read(Field &field_el, std::istream &in_s)
    {
        // Read bigint from stream, reverse bytes in-place and convert to field
        // element.
        if (Form == form_plain) {
            libff::bigint<n> res;
            in_s.read((char *)(&res.data[0]), sizeof(res));
            std::reverse((char *)(&res), (char *)(&res + 1));
            field_el = Field(res);
        } else {
            libff::bigint<n> &res = field_el.mont_repr;
            in_s.read((char *)(&res.data[0]), sizeof(res));
            std::reverse((char *)(&res), (char *)(&res + 1));
        }
    }
};

} // namespace internal

template<typename BigIntT>
void bigint_from_hex(BigIntT &v, const std::string &hex)
{
    hex_to_bytes_reversed(hex, &v.data[0], sizeof(v.data));
}

template<typename BigIntT> std::string bigint_to_hex(const BigIntT &v)
{
    return bytes_to_hex_reversed(&v.data[0], sizeof(v.data));
}

template<encoding_t Enc, form_t Form, typename FieldT>
void field_read(FieldT &v, std::istream &in_s)
{
    internal::field_element_codec<Enc, FieldT>::read(v, in_s);
}

template<encoding_t Enc, form_t Form, typename FieldT>
void field_write(const FieldT &v, std::ostream &out_s)
{
    internal::field_element_codec<Enc, FieldT>::write(v, out_s);
}

} // namespace libff

#endif // __LIBFF_ALGEBRA_FIELDS_FIELD_SERIALIZATION_TCC__
