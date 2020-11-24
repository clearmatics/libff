#!/usr/bin/env sage -python

from sage.all import *
import sys
sys.path.append("../")
import params_generator

# Curve parameter
u = 0x8508c00000000001

def t_0(x):
    return x**5  - 3 *x**4  + 3 *x**3  - x + 3

def y_0(x):
    return t_0(x)/3

def t_1(x):
    return -x**5  + 3 *x**4  - 3 *x**3  + x

def y_1(x):
    return (-t_1(x))/3

def t_i(x, i):
    if i not in {0 , 1 }:
        raise BaseException("i must be 0 or 1")
    if i == 0 :
        return t_0(x)
    return t_1(x)

def y_i(x, i):
    if i not in {0 , 1 }:
        raise BaseException("i must be 0 or 1")
    if i == 0 :
        return y_0(x)
    return y_1(x)

# Lifting cofactors choosen
h_t = 13
h_y = 9

# Choosen set of solutions (i \in {0,1})
i = 0

# See Table 4: https://eprint.iacr.org/2020/351.pdf
def r(x):
    return (x**6  - 2 *x**5  + 2 *x**3  + x + 1 )/3

def t(x):
    return r(x) * h_t + t_i(x, i)

def y(x):
    return r(x) * h_y + y_i(x, i)

def q(x):
    return (t(x)**2  + 3 *y(x)**2 )/4

prime_r = r(u)
#params_generator.generate_libff_Fp_model_params(prime_r)
Fr = GF(prime_r)

prime_q = q(u)
#params_generator.generate_libff_Fp_model_params(prime_q)
Fq = GF(prime_q)


# G1 cofactor
def g1_h(x):
    return ((103*x^6) - (173*x^5) - (96*x^4) + (293*x^3) + (21*x^2) + (52*x) + 172) // 3


def compute_endomorphism_coefficients():
    """
    Return elements \beta of order 3 in Fq, for use in the endomorphism
    \sigma(x, y) = (\beta * x, y).

    See Section 3.1 and 3.2: https://eprint.iacr.org/2020/351.pdf
    """
    Fqx.<x> = PolynomialRing(Fq)
    beta_eqn = x^2 + x + 1
    roots = [root[0] for root in beta_eqn.roots()]

    # Check that \beta_2 corresponds to the "\omega_1" formula from Section 3.1
    v = Fq(u)
    beta_2 = roots[1]
    assert beta_2 == (((103*v^(11)) - (482*v^(10)) + (732*v^9) + (62*v^8) - (1249*v^7) + (1041*v^6) + (214*v^5) - (761*v^4) + (576*v^3) + (11*v^2) - (265*v) + 66) // 21)

    # Return the G1 coefficient "\omega_1" first.
    return roots[1], roots[0]


def g1_fast_subgroup_check_coefficients(x):
    """
    See Section 3.1: https://eprint.iacr.org/2020/351.pdf
    """
    x_2 = x*x
    return x_2, x_2*x


# G2 cofactor
def g2_h(x):
    return ((103*x^6) - (173*x^5) - (96*x^4) + (293*x^3) + (21*x^2) + (52*x) + 151) // 3

print(f"q = {prime_q}")
print(f"r = {prime_r}")

curve = EllipticCurve(Fq, [0, -1])
curve_order = curve.order()
print(f"curve_order = {curve_order}")

m_twist = EllipticCurve(Fq, [0, 4])
twist_order = m_twist.order()
print(f"twist_order = {twist_order}")

h1 = g1_h(u)
print("h1 = {}".format(h1))
h2 = g2_h(u)
print("h2 = {}".format(h2))
print(f"u = {u}")

beta_1, beta_2 = compute_endomorphism_coefficients()
print(f"  \\beta_1 = {beta_1}")
print(f"g1 fast subgroup check:")
print("  [\\sigma_1(x,y)=(\\beta_1 * x, y) on G1]")
print("  [check [u]P+P+\\sigma_1([u_3]P−[u_2]P+P) == 0]")
print(f"g2 fast subgroup check:")
print("  [\\sigma_2(x,y)=(\\beta_1 * x, y) on G2]")
# print(f"  beta_2 = {beta_2}")
print("  [check -[u]P-P+\\sigma_2([u_3]P−[u_2]P-[u]P) == 0]")
