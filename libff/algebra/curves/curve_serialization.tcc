/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef __LIBFF_ALGEBRA_CURVES_CURVE_SERIALIZATION_TCC__
#define __LIBFF_ALGEBRA_CURVES_CURVE_SERIALIZATION_TCC__

#include "libff/algebra/curves/curve_serialization.hpp"

namespace libff
{

template<
    encoding_t Enc,
    form_t Form,
    compression_t Comp,
    typename GroupT>
class group_element_codec;

/// Json encoding/decoding can be performed only with compression_off.
template<form_t Form, typename GroupT>
class group_element_codec<encoding_json, Form, compression_off, GroupT>
{
public:
    static void write(const GroupT &group_el, std::ostream &out_s)
    {
        GroupT affine_p = group_el;
        affine_p.to_affine_coordinates();
        out_s << "[";
        field_write<encoding_json, Form>(affine_p.X, out_s);
        out_s << ",";
        field_write<encoding_json, Form>(affine_p.Y, out_s);
        out_s << "]";
    }
    static void read(GroupT &group_el, std::istream &in_s)
    {
        using base_field = typename std::decay<decltype(group_el.X)>::type;

        char sep;

        in_s >> sep;
        if (sep != '[') {
            throw std::runtime_error(
                "expected opening bracket reading group element");
        }
        field_read<encoding_json, Form>(group_el.X, in_s);

        in_s >> sep;
        if (sep != ',') {
            throw std::runtime_error("expected comma reading group element");
        }

        field_read<encoding_json, Form>(group_el.Y, in_s);
        in_s >> sep;
        if (sep != ']') {
            throw std::runtime_error(
                "expected closing bracket reading group element");
        }

        if (group_el.X.is_zero() && (group_el.Y == base_field::one())) {
            group_el.Z = group_el.Z.zero();
        } else {
            group_el.Z = group_el.Z.one();
        }
    }
};

/// Generic binary uncompressed encoding / decoding
template<
    form_t Form,
    typename GroupT>
class group_element_codec<encoding_binary, Form, compression_off, GroupT>
{
public:
    static void write(const GroupT &group_el, std::ostream &out_s)
    {
        GroupT affine_p = group_el;
        affine_p.to_affine_coordinates();
        field_write<encoding_binary, Form>(affine_p.X, out_s);
        field_write<encoding_binary, Form>(affine_p.Y, out_s);
    }
    static void read(GroupT &group_el, std::istream &in_s)
    {
        using base_field = typename std::decay<decltype(group_el.X)>::type;
        field_read<encoding_binary, Form>(group_el.X, in_s);
        field_read<encoding_binary, Form>(group_el.Y, in_s);
        if (group_el.X.is_zero() && group_el.Y == base_field::one()) {
            group_el.Z = base_field::zero();
        } else {
            group_el.Z = base_field::one();
        }
    }
};

/// Binary compressed encoding / decoding (plain form only). First 2 bits are
/// used to signify GroupT::zero() and -ve Y, respectively.
template<typename GroupT>
class group_element_codec<encoding_binary, form_plain, compression_on, GroupT>
{
public:
    static bool verify_flag_capacity()
    {
        using Fq = typename GroupT::base_field;
        const size_t fq_bits = Fq::num_bits;
        const size_t fq_size_on_disk = 8 * sizeof(Fq);
        return (fq_size_on_disk - fq_bits) >= 2;
    }

    static void write(const GroupT &group_el, std::ostream &out_s)
    {
        assert(verify_flag_capacity());
        using Fq = typename std::decay<decltype(group_el.X)>::type;
        using BigIntT = typename std::decay<decltype(group_el.X.mont_repr)>::type;
        constexpr size_t flag_shift = (8 * sizeof(mp_limb_t)) - 2;
        constexpr mp_limb_t value_mask = (((mp_limb_t)1) << flag_shift) - 1;
        UNUSED(value_mask);

        // Construct X in plain form and write directly
        if (!group_el.is_zero()) {
            GroupT affine(group_el);
            affine.to_affine_coordinates();
            BigIntT v = affine.Y.as_bigint();
            const mp_limb_t flags = v.data[0] & 1;

            // Re-use v to construct the plain X with flags
            v = affine.X.as_bigint();
            assert(
                (v.data[BigIntT::N - 1] & value_mask) == v.data[BigIntT::N - 1]);
            v.data[BigIntT::N - 1] = v.data[BigIntT::N - 1] | (flags << flag_shift);
            std::reverse((char *)(&v), (char *)(&v + 1));
            out_s.write((const char *)(&v.data[0]), sizeof(v));
        } else {
            Fq x;
            x.mont_repr.data[BigIntT::N - 1] = 0x2ull << flag_shift;
            field_write<encoding_binary, form_montgomery>(x, out_s);
        }
    }

    static void read(GroupT &group_el, std::istream &in_s)
    {
        assert(verify_flag_capacity());
        using Fq = typename std::decay<decltype(group_el.X)>::type;
        using BigIntT = typename std::decay<decltype(group_el.X.mont_repr)>::type;
        constexpr size_t flag_shift = (8 * sizeof(mp_limb_t)) - 2;
        constexpr mp_limb_t value_mask = (((mp_limb_t)1) << flag_shift) - 1;

        // Reading the value raw, and extract the flags
        Fq X_and_flags;
        field_read<encoding_binary, form_montgomery>(X_and_flags, in_s);
        const mp_limb_t hi_limb = X_and_flags.mont_repr.data[BigIntT::N - 1];
        const mp_limb_t flags = hi_limb >> flag_shift;

        if (0 == (flags & 0x2)) {
            // Remove the flags from X_and_flags, and use to create X (in
            // Montgomery form).
            X_and_flags.mont_repr.data[BigIntT::N - 1] = hi_limb & value_mask;
            group_el.X = Fq(X_and_flags.mont_repr);
            group_el.Y = curve_point_y_at_x<GroupT>(group_el.X);

            // Reuse X_and_flags.mont_repr to hold the reduced Y value. Invert
            // Y if least-significant bit does not match the flag.
            X_and_flags.mont_repr = group_el.Y.as_bigint();
            if ((flags & 1) != (X_and_flags.mont_repr.data[0] & 1)) {
                group_el.Y = - group_el.Y;
            }

            group_el.Z = Fq::one();
        } else {
            group_el = GroupT::zero();

        }
    }
};

template<
    encoding_t Enc,
    form_t Form,
    compression_t Comp,
    typename GroupT>
void group_read(GroupT &v, std::istream &in_s)
{
    group_element_codec<Enc, Form, Comp, GroupT>::read(v, in_s);
}

template<
    encoding_t Enc,
    form_t Form,
    compression_t Comp,
    typename GroupT>
void group_write(const GroupT &v, std::ostream &out_s)
{
    group_element_codec<Enc, Form, Comp, GroupT>::write(v, out_s);
}

} // namespace libff

#endif // __LIBFF_ALGEBRA_CURVES_CURVE_SERIALIZATION_TCC__
