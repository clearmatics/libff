/** @file
 *****************************************************************************
 Declaration of arithmetic in the finite field F[(p^2)^3]
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP6_3OVER2_HPP_
#define FP6_3OVER2_HPP_
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp2.hpp>
#include <vector>

namespace libff
{

template<mp_size_t n, const bigint<n> &modulus> class Fp6_3over2_model;

template<mp_size_t n, const bigint<n> &modulus>
std::ostream &operator<<(std::ostream &, const Fp6_3over2_model<n, modulus> &);

template<mp_size_t n, const bigint<n> &modulus>
std::istream &operator>>(std::istream &, Fp6_3over2_model<n, modulus> &);

/// Arithmetic in the finite field F[(p^2)^3].
///
/// Let p := modulus. This interface provides arithmetic for the extension
/// field Fp6 = Fp2[V]/(V^3-non_residue) where non_residue is in Fp.
///
/// ASSUMPTION: p = 1 (mod 6)
template<mp_size_t n, const bigint<n> &modulus> class Fp6_3over2_model
{
public:
    typedef Fp_model<n, modulus> my_Fp;
    typedef Fp2_model<n, modulus> my_Fp2;

    // (modulus^6-1)/2
    static bigint<6 * n> euler;
    // modulus^6 = 2^s * t + 1
    static std::size_t s;
    // with t odd
    static bigint<6 * n> t;
    // (t-1)/2
    static bigint<6 * n> t_minus_1_over_2;
    // a quadratic nonresidue in Fp6
    static Fp6_3over2_model<n, modulus> nqr;
    // nqr^t
    static Fp6_3over2_model<n, modulus> nqr_to_t;

    static my_Fp2 non_residue;
    /// non_residue^((modulus^i-1)/3)   for i=0,1,2,3,4,5
    static my_Fp2 Frobenius_coeffs_c1[6];
    /// non_residue^((2*modulus^i-2)/3) for i=0,1,2,3,4,5
    static my_Fp2 Frobenius_coeffs_c2[6];

    static const size_t tower_extension_degree = 3;

    my_Fp2 coeffs[3];
    Fp6_3over2_model(){};
    Fp6_3over2_model(const my_Fp2 &c0, const my_Fp2 &c1, const my_Fp2 &c2)
    {
        this->coeffs[0] = c0;
        this->coeffs[1] = c1;
        this->coeffs[2] = c2;
        return;
    };

    void clear()
    {
        coeffs[0].clear();
        coeffs[1].clear();
        coeffs[2].clear();
    }
    void print() const
    {
        printf("c0/c1/c2:\n");
        coeffs[0].print();
        coeffs[1].print();
        coeffs[2].print();
    }

    static Fp6_3over2_model<n, modulus> zero();
    static Fp6_3over2_model<n, modulus> one();
    static Fp6_3over2_model<n, modulus> random_element();

    bool is_zero() const
    {
        return coeffs[0].is_zero() && coeffs[1].is_zero() &&
               coeffs[2].is_zero();
    }
    bool operator==(const Fp6_3over2_model &other) const;
    bool operator!=(const Fp6_3over2_model &other) const;

    Fp6_3over2_model operator+(const Fp6_3over2_model &other) const;
    Fp6_3over2_model operator-(const Fp6_3over2_model &other) const;
    Fp6_3over2_model operator*(const Fp6_3over2_model &other) const;
    Fp6_3over2_model operator-() const;
    Fp6_3over2_model squared() const;
    Fp6_3over2_model inverse() const;
    Fp6_3over2_model Frobenius_map(unsigned long power) const;

    static my_Fp2 mul_by_non_residue(const my_Fp2 &elt);

    template<mp_size_t m>
    Fp6_3over2_model operator^(const bigint<m> &other) const;

    static bigint<n> base_field_char() { return modulus; }
    static constexpr size_t extension_degree() { return 6; }

    friend std::ostream &operator<<<n, modulus>(
        std::ostream &out, const Fp6_3over2_model<n, modulus> &el);
    friend std::istream &operator>>
        <n, modulus>(std::istream &in, Fp6_3over2_model<n, modulus> &el);
};

template<mp_size_t n, const bigint<n> &modulus>
std::ostream &operator<<(
    std::ostream &out, const std::vector<Fp6_3over2_model<n, modulus>> &v);

template<mp_size_t n, const bigint<n> &modulus>
std::istream &operator>>(
    std::istream &in, std::vector<Fp6_3over2_model<n, modulus>> &v);

template<mp_size_t n, const bigint<n> &modulus>
Fp6_3over2_model<n, modulus> operator*(
    const Fp_model<n, modulus> &lhs, const Fp6_3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n> &modulus>
Fp6_3over2_model<n, modulus> operator*(
    const Fp2_model<n, modulus> &lhs, const Fp6_3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n> &modulus>
bigint<6 * n> Fp6_3over2_model<n, modulus>::euler;

template<mp_size_t n, const bigint<n> &modulus>
size_t Fp6_3over2_model<n, modulus>::s;

template<mp_size_t n, const bigint<n> &modulus>
bigint<6 * n> Fp6_3over2_model<n, modulus>::t;

template<mp_size_t n, const bigint<n> &modulus>
bigint<6 * n> Fp6_3over2_model<n, modulus>::t_minus_1_over_2;

template<mp_size_t n, const bigint<n> &modulus>
Fp6_3over2_model<n, modulus> Fp6_3over2_model<n, modulus>::nqr;

template<mp_size_t n, const bigint<n> &modulus>
Fp6_3over2_model<n, modulus> Fp6_3over2_model<n, modulus>::nqr_to_t;

template<mp_size_t n, const bigint<n> &modulus>
Fp2_model<n, modulus> Fp6_3over2_model<n, modulus>::non_residue;

template<mp_size_t n, const bigint<n> &modulus>
Fp2_model<n, modulus> Fp6_3over2_model<n, modulus>::Frobenius_coeffs_c1[6];

template<mp_size_t n, const bigint<n> &modulus>
Fp2_model<n, modulus> Fp6_3over2_model<n, modulus>::Frobenius_coeffs_c2[6];

} // namespace libff

#include <libff/algebra/fields/fp6_3over2.tcc>

#endif // FP6_3OVER2_HPP_
