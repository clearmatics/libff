/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef CURVE_UTILS_HPP_
#define CURVE_UTILS_HPP_
#include <cstdint>
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/bigint.hpp>

namespace libff
{

// The field element (base or extension field) type for the coordinates of any
// group.
template<typename GroupT>
using group_coord_type =
    typename std::decay<decltype(((GroupT *)nullptr)->X)>::type;

template<typename GroupT, mp_size_t m>
GroupT scalar_mul(const GroupT &base, const bigint<m> &scalar);

// Utility function to compute Y coordinate of a point on the curve E(Fq) with
// the given x coordinate. This function does not check whether E(Fq) has a
// solution at x, and will hang indefinitely if it does not.
template<typename GroupT>
group_coord_type<GroupT> curve_point_y_at_x(const group_coord_type<GroupT> &x);

// Utility function to compute a point on the curve E(Fq) with the given x
// coordinate. If the curve has no solution at the given x, this function
// throws an exception.
template<typename GroupT>
GroupT g1_curve_point_at_x(const typename GroupT::base_field &x);

// Utility function to compute a point on the twisted curve E'(Fqe) with the
// given x coordinate.
template<typename GroupT>
GroupT g2_curve_point_at_x(const typename GroupT::twist_field &x);

} // namespace libff
#include <libff/algebra/curves/curve_utils.tcc>

#endif // CURVE_UTILS_HPP_
