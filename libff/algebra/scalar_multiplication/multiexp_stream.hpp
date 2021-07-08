/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MULTIEXP_STREAM_HPP_
#define MULTIEXP_STREAM_HPP_

#include "libff/algebra/serialization.hpp"

#include <cstddef>
#include <vector>

namespace libff
{

/// Read base elements from a stream. More intermediate memory is used (offset
/// by the fact that base elements are streamed and therefore not all
/// memory-resident) to reduce the number of internal passes. Currently
/// processing is single-threaded (although element streaming happens in a
/// separate temporary thread).
template<form_t Form, compression_t Comp, typename GroupT, typename FieldT>
GroupT multi_exp_stream(
    std::istream &base_elements_in, const std::vector<FieldT> &exponents);

} // namespace libff

#include "libff/algebra/scalar_multiplication/multiexp_stream.tcc"

#endif // MULTIEXP_STREAM_HPP_
