#!/usr/bin/env sage -python

from sage.all import *
import sys
sys.path.append("../")
import params_generator

# Curve parameter
u = 0x8508c00000000001

# Prime order of the subgroup we work in
def r(x):
    return (x^4) - (x^2) + 1

# Prime used to generate the base finite field
def q(x):
    return ((x - 1)^2) * (r(x) / 3) + x

prime_r = r(u)
assert(prime_r == 8444461749428370424248824938781546531375899335154063827935233455917409239041)
prime_q = q(u)
assert(prime_q == 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177)
if (mod(prime_q, 6) != 1):
    raise BaseException("Unexpected: q should be = 1 (mod 6). See: https://eprint.iacr.org/2007/390.pdf")

print('prime_r = {}'.format(prime_r))
params_generator.generate_libff_Fp_model_params(prime_r)
Fr = GF(prime_r)

print('prime_q = {}'.format(prime_q))
params_generator.generate_libff_Fp_model_params(prime_q)
Fq = GF(prime_q)

# Generate Fp2_model
euler_fp2 = (prime_q**2 - 1)/2; print('euler = {}'.format(euler_fp2))

# factorization = factor(prime_q**2 - 1)
# t = 0
# term_2 = factorization[0]
# counter = 0
# if term_2[0] != 2:
#     raise BaseException("The prime decomposition doesn't have any factor 2 - the 'high 2-adicity' requirement isn't respected")
# while not(is_odd(t)):
#     s = term_2[1] - counter;
#     t = (prime_q**2 - 1)/(2**s);
#     counter = counter + 1
#
# print('s = {}'.format(s))
# is_odd(t); print('t = {}'.format(t))

# The t and s parameters below are the result of the code commented above. It takes ages to run so the code has been commented out
# and the result given directly to avoid running expensive code all the time
t = 475404855284145089315325463221726483993816145966867441829193658311651761271425728823393990805904040047516478740222806302278755994777496288961383541476974255391881599499962735436887347234371823579436839914935817251
s = 47

t_minus_1_over_2 = (t-1)/2
print('t_minus_1_over_2 = {}'.format(t_minus_1_over_2))

# Fq2 is constructed as Fq[u]/(u^2 - beta)Fq[u], where beta = -5
non_residue = Fq(-5)
print('non_residue = {}'.format(non_residue))
Fqx.<j> = PolynomialRing(Fq, 'j')
assert(Fqx(j^2 + 5).is_irreducible())
Fq2.<u> = GF(prime_q^2, modulus=j^2 + 5)

nqr = Fq2(u)
assert(nqr^euler_fp2 == Fq2(-1))
nqr_to_t = pow(nqr, t, Fq2(u^2 + 5))
print('nqr_to_t = {}'.format(nqr_to_t))

m_twist = EllipticCurve(Fq2, [0, (1 * u)])
d_twist = EllipticCurve(Fq2, [0, (1 / u)])
if (m_twist.order() % prime_r) == 0:
    print "Good twist is: M-type"
elif (d_twist.order() % prime_r) == 0:
    print "Good twist is: D-type"
else:
    raise BaseException("Error. None of the proposed twists has order divisible by r")

frob_coeff_c1_1 = non_residue^((0/2)) # = 1
frob_coeff_c1_2 = non_residue^(((prime_q) - 1)/2)

### Fq6 is constructed as Fq2[v]/(v^3 - u)Fq2[v]
##Fq2x.<l> = PolynomialRing(Fq2, 'l')
##assert(Fq2x(l^3 - u).is_irreducible())
###Fq6.<v> = GF((prime_q^2)^3, modulus=l^3 - u)
##non_residue = Fq6(v);
##print('non_residue = {}'.format(non_residue))
##
### Fq12 is constructed as Fq6[w]/(w^2 + v)F[w]
##Fq3x.<w> = PolynomialRing(Fq3, 'w')
##assert(Fqx(w^2 + v).is_irreducible())
##Fq12.<n> = GF(((prime_q^2)^3)^2, modulus=w^2 +