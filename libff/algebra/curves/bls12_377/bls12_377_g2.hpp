/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BLS12_377_G2_HPP_
#define BLS12_377_G2_HPP_
#include <vector>

#include <libff/algebra/curves/bls12_377/bls12_377_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class bls12_377_G2;
std::ostream& operator<<(std::ostream &, const bls12_377_G2&);
std::istream& operator>>(std::istream &, bls12_377_G2&);

class bls12_377_G2 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static bls12_377_G2 G2_zero;
    static bls12_377_G2 G2_one;
    static bls12_377_Fq2 coeff_a;
    static bls12_377_Fq2 coeff_b;

    typedef bls12_377_Fq base_field;
    typedef bls12_377_Fq2 twist_field;
    typedef bls12_377_Fr scalar_field;

    // Cofactor
    static const mp_size_t h_bitcount = 502;
    static const mp_size_t h_limbs = (h_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
    static bigint<h_limbs> h;

    bls12_377_Fq2 X, Y, Z;

    // using Jacobian coordinates
    bls12_377_G2();
    bls12_377_G2(const bls12_377_Fq2& X, const bls12_377_Fq2& Y, const bls12_377_Fq2& Z) : X(X), Y(Y), Z(Z) {};

    static bls12_377_Fq2 mul_by_b(const bls12_377_Fq2 &elt);

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const bls12_377_G2 &other) const;
    bool operator!=(const bls12_377_G2 &other) const;

    bls12_377_G2 operator+(const bls12_377_G2 &other) const;
    bls12_377_G2 operator-() const;
    bls12_377_G2 operator-(const bls12_377_G2 &other) const;

    bls12_377_G2 add(const bls12_377_G2 &other) const;
    bls12_377_G2 mixed_add(const bls12_377_G2 &other) const;
    bls12_377_G2 dbl() const;
    bls12_377_G2 mul_by_q() const;
    bls12_377_G2 mul_by_cofactor() const;

    bool is_well_formed() const;
    bool is_in_safe_subgroup() const;

    static bls12_377_G2 zero();
    static bls12_377_G2 one();
    static bls12_377_G2 random_element();

    static size_t size_in_bits() { return twist_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    void write_uncompressed(std::ostream &) const;
    void write_compressed(std::ostream &) const;
    static void read_uncompressed(std::istream &, bls12_377_G2 &);
    static void read_compressed(std::istream &, bls12_377_G2 &);
    static void read_uncompressed_unsafe(std::istream &, bls12_377_G2 &);
    static void read_compressed_unsafe(std::istream &, bls12_377_G2 &);

    static void batch_to_special_all_non_zeros(std::vector<bls12_377_G2> &vec);
};

template<mp_size_t m>
bls12_377_G2 operator*(const bigint<m> &lhs, const bls12_377_G2 &rhs)
{
    return scalar_mul<bls12_377_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
bls12_377_G2 operator*(const Fp_model<m,modulus_p> &lhs, const bls12_377_G2 &rhs)
{
    return scalar_mul<bls12_377_G2, m>(rhs, lhs.as_bigint());
}


} // libff
#endif // BLS12_377_G2_HPP_
