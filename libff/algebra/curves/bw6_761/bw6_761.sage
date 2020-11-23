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

def g1_beta(x):
    """
    G1 \beta (element order 3 in Fq) for endomorphism \sigma(x, y) = (\beta * x, y).
    See Section 3.1: https://eprint.iacr.org/2020/351.pdf
    """
    x = Fq(x)
    return (((103*x^(11)) - (482*x^(10)) + (732*x^9) + (62*x^8) - (1249*x^7) + (1041*x^6) + (214*x^5) - (761*x^4) + (576*x^3) + (11*x^2) - (265*x) + 66) // 21)

def g1_fast_subgroup_check_coefficients(x):
    """
    See Section 3.1: https://eprint.iacr.org/2020/351.pdf
    """
    x_2 = x*x
    return x_2, x_2*x

# G2 cofactor
def g2_h(x):
    return ((103*x^6) - (173*x^5) - (96*x^4) + (293*x^3) + (21*x^2) + (52*x) + 151) // 3

h1 = g1_h(u)
print('h1 = {}'.format(h1))
h2 = g2_h(u)
print('h2 = {}'.format(h2))

print(f"g1 fast subgroup check:")
beta = g1_beta(u)
print("  [\sigma(x,y)=(\beta * x, y)]")
print(f"  beta = {beta} ")
u_2, u_3 = g1_fast_subgroup_check_coefficients(u)
print("  [check [u]P+P+\sigma([u_3]P−[u_2]P+P) == 0]")
print(f"  u={u},u_2={u_2},u_3={u_3}")
