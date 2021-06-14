/** @file
 *
 * This is a modified version of code that originaly appeared in
 * https://github.com/clearmatics/zeth.
 *
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

namespace internal
{

// Generic class to implement group element read and write methods for various
// combinations of parameters. Expected to define at least methods:
//
//   static void write(const GroupT &group_el, std::ostream &out_s);
//   static void read(GroupT &group_el, std::istream &in_s);
template<encoding_t Enc, form_t Form, compression_t Comp, typename GroupT>
class group_element_codec;

// Json encoding/decoding. Supported only for compression_off.
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

// Generic uncompressed binary encoding / decoding
template<form_t Form, typename GroupT>
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

// Binary compressed encoding / decoding. The base field must allow for 2 flags
// bits to fit in the high-order limb:
//   bit 0 (mask 0x1) : Lowest order bit of Y coord
//   bit 1 (mask 0x2) : First 2 bits signify GroupT::zero()
//
// Note, we rely on the field_write_with_flags call to ensure capacity in the
// limb and to shift the flags to the appropriate place.
template<form_t Form, typename GroupT>
class group_element_codec<encoding_binary, Form, compression_on, GroupT>
{
public:
    static void write(const GroupT &group_el, std::ostream &out_s)
    {
        using BigIntT =
            typename std::decay<decltype(group_el.X.mont_repr)>::type;

        if (!group_el.is_zero()) {
            GroupT affine(group_el);
            affine.to_affine_coordinates();

            BigIntT v = affine.Y.as_bigint();
            const mp_limb_t flags = v.data[0] & 1;

            field_write_with_flags<encoding_binary, Form>(
                affine.X, flags, out_s);
        } else {
            // Use Montgomery encoding, to avoid wasting time reducing.
            field_write_with_flags<encoding_binary, form_montgomery>(
                group_el.X, 0x2, out_s);
        }
    }
    static void read(GroupT &group_el, std::istream &in_s)
    {
        using Fq = typename std::decay<decltype(group_el.X)>::type;

        mp_limb_t flags;
        Fq X;

        // Initially, read the value without reduction. Only reduce if
        // necessary. (For the Form == form_montgomery case, we could read
        // directly into group_el.X. This may be worth measuring.)

        field_read_with_flags<encoding_binary, form_montgomery>(X, flags, in_s);
        if (0 == (flags & 0x2)) {
            if (Form == form_plain) {
                group_el.X = Fq(X.mont_repr);
            } else {
                group_el.X = X;
            }
            group_el.Y = curve_point_y_at_x<GroupT>(group_el.X);

            // Reuse X.mont_repr to hold the reduced Y value. Invert
            // Y if least-significant bit does not match the flag.
            X.mont_repr = group_el.Y.as_bigint();
            if ((flags & 1) != (X.mont_repr.data[0] & 1)) {
                group_el.Y = -group_el.Y;
            }

            group_el.Z = Fq::one();
        } else {
            group_el = GroupT::zero();
        }
    }
};

} // namespace internal

template<encoding_t Enc, form_t Form, compression_t Comp, typename GroupT>
void group_read(GroupT &v, std::istream &in_s)
{
    internal::group_element_codec<Enc, Form, Comp, GroupT>::read(v, in_s);
}

template<encoding_t Enc, form_t Form, compression_t Comp, typename GroupT>
void group_write(const GroupT &v, std::ostream &out_s)
{
    internal::group_element_codec<Enc, Form, Comp, GroupT>::write(v, out_s);
}

} // namespace libff

#endif // __LIBFF_ALGEBRA_CURVES_CURVE_SERIALIZATION_TCC__