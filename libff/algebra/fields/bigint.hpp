/** @file
 *****************************************************************************
 Declaration of bigint wrapper class around GMP's MPZ long integers.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BIGINT_HPP_
#define BIGINT_HPP_
#include <cstddef>
#include <gmp.h>
#include <iostream>
#include <libff/common/serialization.hpp>

namespace libff
{

template<mp_size_t n> class bigint;
template<mp_size_t n>
std::ostream &operator<<(std::ostream &, const bigint<n> &);
template<mp_size_t n> std::istream &operator>>(std::istream &, bigint<n> &);

/// Wrapper class around GMP's MPZ long integers. It supports arithmetic
/// operations, serialization and randomization. Serialization is fragile, see
/// common/serialization.hpp.
template<mp_size_t n> class bigint
{
public:
    static const mp_size_t N = n;

    mp_limb_t data[n] = {0};

    bigint() = default;
    /// Initalize from a small integer
    bigint(const unsigned long x);
    /// Initialize from a string containing an integer in decimal notation
    bigint(const char *s);
    /// Initialize from MPZ element
    bigint(const mpz_t r);

    void print() const;
    void print_hex() const;
    bool operator==(const bigint<n> &other) const;
    bool operator!=(const bigint<n> &other) const;
    void clear();
    bool is_zero() const;

    /// The number of bits representable by this bigint type
    static constexpr size_t max_bits() { return n * GMP_NUMB_BITS; };

    /// The number of bits in this specific bigint value, i.e., position of the
    /// most-significant 1
    size_t num_bits() const;

    unsigned long as_ulong() const; /// Return the last limb of the integer
    void to_mpz(mpz_t r) const;
    bool test_bit(const std::size_t bitno) const;

    bigint &randomize();

    friend std::ostream &operator<<<n>(std::ostream &out, const bigint<n> &b);
    friend std::istream &operator>><n>(std::istream &in, bigint<n> &b);
};

} // namespace libff

#include <libff/algebra/fields/bigint.tcc>

#endif // BIGINT_HPP_
