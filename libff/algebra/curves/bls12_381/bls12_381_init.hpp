/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BLS12_381_INIT_HPP_
#define BLS12_381_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
//#include <libff/algebra/curves/bls12_381/bls12_381_fields.hpp> (VV)
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp12_2over3over2.hpp>
#include <libff/algebra/fields/fp2.hpp>
#include <libff/algebra/fields/fp6_3over2.hpp>

namespace libff
{

const mp_size_t bls12_381_r_bitcount = 255;
const mp_size_t bls12_381_q_bitcount = 381;

const mp_size_t bls12_381_r_limbs =
    (bls12_381_r_bitcount + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
const mp_size_t bls12_381_q_limbs =
    (bls12_381_q_bitcount + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;

extern bigint<bls12_381_q_limbs> bls12_381_modulus_q;
// by example of bls12_377_modulus_q in bls12_377_init.hpp
// extern bigint<bls12_381_q_limbs> bw6_761_modulus_r; // (VV)

extern bigint<bls12_381_r_limbs> bls12_381_modulus_r;

// by example of bls12_377_modulus_q in bls12_377_init.hpp
//#define bls12_381_modulus_q bw6_761_modulus_r // (VV)

typedef Fp_model<bls12_381_r_limbs, bls12_381_modulus_r> bls12_381_Fr;
typedef Fp_model<bls12_381_q_limbs, bls12_381_modulus_q> bls12_381_Fq;
typedef Fp2_model<bls12_381_q_limbs, bls12_381_modulus_q> bls12_381_Fq2;
typedef Fp6_3over2_model<bls12_381_q_limbs, bls12_381_modulus_q> bls12_381_Fq6;
typedef Fp12_2over3over2_model<bls12_381_q_limbs, bls12_381_modulus_q>
    bls12_381_Fq12;
typedef bls12_381_Fq12 bls12_381_GT;

// parameters for the curve E/Fq : y^2 = x^3 + b
extern bls12_381_Fq bls12_381_coeff_b;
extern bigint<bls12_381_r_limbs>
    bls12_381_trace_of_frobenius; // from bls12_381 (VV)
// parameters for the twisted curve E'/Fq2 : y^2 = x^3 + b/xi
extern bls12_381_Fq2 bls12_381_twist;
extern bls12_381_Fq2 bls12_381_twist_coeff_b;
extern bls12_381_Fq bls12_381_twist_mul_by_b_c0;
extern bls12_381_Fq bls12_381_twist_mul_by_b_c1;
extern bls12_381_Fq2 bls12_381_twist_mul_by_q_X;
extern bls12_381_Fq2 bls12_381_twist_mul_by_q_Y;

// Coefficient \beta in endomorphism (x, y) -> (\beta * x, y) (from bls12_377
// VV)
extern bls12_381_Fq bls12_381_g1_endomorphism_beta;
extern bigint<bls12_381_r_limbs> bls12_381_g1_safe_subgroup_check_c1;

// Coefficients for G2 untwist-frobenius-twist (from bls12_377 VV)
extern bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_v;
extern bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_w_3;
extern bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_v_inverse;
extern bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_w_3_inverse;

// parameters for pairing
extern bigint<bls12_381_q_limbs> bls12_381_ate_loop_count;
extern bool bls12_381_ate_is_loop_count_neg;
extern bigint<12 * bls12_381_q_limbs> bls12_381_final_exponent;
extern bigint<bls12_381_q_limbs> bls12_381_final_exponent_z;
extern bool bls12_381_final_exponent_is_z_neg;

void init_bls12_381_params();

// void init_bls12_381_fields(); // from scipr-lab (VV)

class bls12_381_G1;
class bls12_381_G2;

} // namespace libff
#endif // BLS12_381_INIT_HPP_
