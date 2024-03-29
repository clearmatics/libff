/** @file
*****************************************************************************
* @author     This file is part of libff, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef BLS12_377_PP_HPP_
#define BLS12_377_PP_HPP_
#include <libff/algebra/curves/bls12_377/bls12_377_g1.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_g2.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_init.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff
{

class bls12_377_pp
{
public:
    static const std::string name;

    typedef bls12_377_Fr Fp_type;
    typedef bls12_377_G1 G1_type;
    typedef bls12_377_G2 G2_type;
    typedef bls12_377_G1_precomp G1_precomp_type;
    typedef bls12_377_G2_precomp G2_precomp_type;
    typedef bls12_377_Fq Fq_type;
    typedef bls12_377_Fq2 Fqe_type;
    typedef bls12_377_Fq12 Fqk_type;
    typedef bls12_377_GT GT_type;

    static const bool has_affine_pairing = false;

    static void init_public_params();
    static bls12_377_GT final_exponentiation(const bls12_377_Fq12 &elt);
    static bls12_377_G1_precomp precompute_G1(const bls12_377_G1 &P);
    static bls12_377_G2_precomp precompute_G2(const bls12_377_G2 &Q);
    static bls12_377_Fq12 miller_loop(
        const bls12_377_G1_precomp &prec_P, const bls12_377_G2_precomp &prec_Q);
    static bls12_377_Fq12 double_miller_loop(
        const bls12_377_G1_precomp &prec_P1,
        const bls12_377_G2_precomp &prec_Q1,
        const bls12_377_G1_precomp &prec_P2,
        const bls12_377_G2_precomp &prec_Q2);
    static bls12_377_Fq12 pairing(const bls12_377_G1 &P, const bls12_377_G2 &Q);
    static bls12_377_Fq12 reduced_pairing(
        const bls12_377_G1 &P, const bls12_377_G2 &Q);
};

} // namespace libff

#endif // BLS12_377_PP_HPP_
