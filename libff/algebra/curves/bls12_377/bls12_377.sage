#!/usr/bin/env sage -python3

from sage.all import *
import sys
sys.path.append("../")
import params_generator

# Prime order of the subgroup we work in
def r(x):
    return (x^4) - (x^2) + 1

# Prime used to generate the base finite field
def q(x):
    return ((x - 1)^2) * (r(x) / 3) + x

# Compute G1 cofactor
def g1_h(x):
    return ((x - 1)^2) // 3

def g1_t(x):
    return (x + 1)

# Compute G2 cofactor
# See: Proposition 2, Section 3.4: https://eprint.iacr.org/2015/247.pdf
def g2_h(x):
    return ((x**8) - (4 * (x**7)) + (5 * (x**6)) - (4 * (x**4)) + (6 * (x**3)) - (4 * (x**2)) - (4*x) + 13) // 9

# Computes the order of G1, the safe subgroup of E/Fq
def g1_order(curve_order):
    decomposition = factor(curve_order)
    # Factor returns the prime decomposition and orders prime
    # factors from smaller to biggest
    biggest_factor = decomposition[-1]
    assert(biggest_factor[1] == 1)
    return biggest_factor[0]


def g1_fast_subgroup_check_coefficients(g1_order, q, r):
    """
    Find \beta Fq with multiplicative order 3. For each element (x,y) in G1
    with n = order or (x, y), the endomorphism \sigma(x, y) -> (\beta x, y) is
    equivalent to [\lambda_n], where \lambda_n is a root of x^2 + x + 1 in Fn.

    Find c0 and c1 s.t. c0 + c1 * \lambda_r = r
    (i.e. [c0]P + [c1]\sigma(P) = [r]P for P in r-torsion)

    Ensure that c0 + c1 * x has no roots mod any other subgroup order
    (i.e. [c0]P + [c1]\sigma(P) does not kill members of any other torsion).

    In this case, [c0]P + [c1]\sigma(P) == 0 is a fast r-torsion membership
    check for elements of E(Fq).
    """

    print(f" [r1_0]P + [r2_1]\sigma(P) == 0")

    def _find_lambda_n(Fn):
        """
        Return solutions of \lambda^2 + \lambda + 1 = 0 mod n
        """
        Fny.<y> = PolynomialRing(Fn)
        lambda_r_poly = y^2 + y + 1
        lambda_r_poly_roots = lambda_r_poly.roots()
        return [root[0] for root in lambda_r_poly_roots]

    def _check_lambdas_not_root_in_subgroup_size(n, c0, c1):
        """
        Check that no possible values of lambda_n are roots of c0 + c1*x in GF(n)
        """
        Fn = GF(n)
        lambdas = _find_lambda_n(Fn)
        for l in lambdas:
            if 0 == Fn(c0 + c1 * l):
                raise Exception(f"elements order {n} may be mapped to 0")

    def _check_lambdas_not_root_in_subgroups(group_order, q, r, c0, c1):
        """
        For each prime-power factor n of group_order, perform the check that
        [c0]P + [c1]sigma(P) does not kill the n-torsion.
        """

        def _no_element_order_4():
            """
            E(Fq) contains no elements of order 4 (and hence of order 2^m, m > 1).
              E: y^2 - x^3 - 1
              Any order 2 element has dE/dy == 0, i.e. y == 0.
              Any order 4 element P must thereby have:
                (P + P) = R
                Ry == 0
              By group laws:
                21 * Px^6 - 5 * Px^3 + 2 == 0
            """
            Fq = GF(q)
            Fqp.<x> = PolynomialRing(Fq)
            f = 21 * x^6 - 5 * x^3 + 2
            if len(f.roots()):
                raise Exception(f"E(Fq) may contain elements of order 4")
        _no_element_order_4()

        # For other prime power factors n, show that no lambda_n cannot be a root of:
        #   c0 + c1*x mod n
        for (prime, power) in factor(group_order):
            # print(f"  group_order has factor (p={prime}, pow={power})")
            if prime == r:
                # print("   skipping")
                continue
            if prime == 2:
                # No need to check powers greater than 1
                power = 1

            for p in range(1, power+1):
                n = prime^p
                # print(f"   p={p}, n={n}")
                _check_lambdas_not_root_in_subgroup_size(n, c0, c1)

    # Find an element of Fq with multiplicative order 3
    beta = min(_find_lambda_n(GF(q)))
    print(f" beta: {beta}  [in sigma(x,y) = (beta * x, y)]")

    # lambda for the r-torsion
    lambda_r_poly_roots = _find_lambda_n(GF(r));
    # print(f" lambda_r_poly_roots: {lambda_r_poly_roots}")
    lambda_r = min(lambda_r_poly_roots) # [root[0] for root in lambda_r_poly_roots])
    print(f" (lambda_r: {lambda_r})")

    # Find r1_0, r1_1 s.t. r = r1_0 + r1_1*lambda_r
    r1_1 = int(int(r) // int(lambda_r))
    r1_0 = int(r % r1_1)
    _check_lambdas_not_root_in_subgroups(g1_order, q, r, r1_0, r1_1)

    return r1_0, r1_1


def fuentes_g2_fast_cofactor_coefficients(x):
    """
    See Section 4.1: https://eprint.iacr.org/2017/419.pdf
      [m * h2]P = [x^2 - x - 1]P + [x - 1]ψ(P) + ψ^2(2P)
    where:
      m = 3x^3 - 3
    """
    x_2 = x * x
    c0 = x_2 - x - 1
    c1 = x - 1
    print(f" m = {3*(x_2*x - 1)}")
    return c0, c1


def g2_fast_cofactor_coefficients(curve_order, q, h2):
    """
    Compute faster cofactor multiplication coefficients:
    See Section 3.1: https://eprint.iacr.org/2017/419.pdf

    Multiplication by 3*h2 of P in E'/F_{q^2} can be reduced to:
      [x^3−x^2−x+4]P+ [x^3−x^2−x+1]ψ(P) + [−x^2+2x−1]ψ^2(P)
    """

    # For P in E'(Fq2),
    #   [q]P = [t]ψ(P) - ψ^2(P)
    # We seek
    #   h2_0 and h2_1 s.t. h2 = h2_0 + h2_1 * q mod #E'(Fp2)
    # so that:
    #   [h]P = [h_0]P + [h_1]([t]ψ(P) - ψ^2(P))
    print(f" [h]P = [h2_0]P + [h2_1]([t]ψ(P) - ψ^2(P))")

    # Section 3.1: https://eprint.iacr.org/2017/419.pdf:
    #   x_2 = x*x
    #   x_3 = x * x_2
    #   c0 = x_3 - x_2 - x + 4
    #   c1 = x_3 - x_2 - x + 1
    #   c2 = 2*x - x_2 - 1
    #   return c0, c1, c2
    # This trivial calculation reveals the same coefficients / 3 (resulting
    # in slightly fewer operations)

    h2_1 = int(floor(h2 / q))
    h2_0 = h2 - (h2_1 * q)
    return h2_0, h2_1


def g2_fast_subgroup_check_coefficients(g2_order, q, r, t1, h2):
    """
    Compute coefficients of P, (ψ(P) - P) and ψ^2(P) to quickly multiply by
    some factor k of r, s.t. gcd(k, curve_order/r) == 1. [kr] will kill the
    r-torsion without killing any other components.
    """
    # For P in E'(Fq2), ψ satisfies:
    #   [q1]P = [t1]ψ(P) - ψ^2(P)
    # where ψ is the untwist-frobenius-twist endomorphism

    # Trivial solution:
    #   h1.r = n1 = q1 - t1 + 1
    #   [h1.r] P = [t1]ψ(P) - ψ^2(P) - [t]P + P
    #            = P + [t](ψ(P) - P) - ψ^2(P)
    print(f" 0 == [r2_0]P + [r2_1](ψ(P) - P) + [r2_2]ψ^2(P)")
    g1_order = q - t1 + 1
    h1 = g1_order / r

    # Assert that [r] will only kill the G2 component.
    assert h2 % r != 0
    # Assert that [h1] will not kill any element of G2.
    assert 1 == gcd(h1, h2)

    return 1, t1, -1


def main():
    print("Generating parameters for BLS12-377")
    # Curve parameter
    param = 0x8508c00000000001

    prime_r = r(param)
    assert(prime_r == 8444461749428370424248824938781546531375899335154063827935233455917409239041)

    prime_q = q(param)
    assert(prime_q == 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177)
    if (mod(prime_q, 6) != 1):
        raise BaseException("Unexpected: q should be = 1 (mod 6). See: https://eprint.iacr.org/2007/390.pdf")

    # Scalar field
    print('prime_r = {}'.format(prime_r))
    #params_generator.generate_libff_Fp_model_params(prime_r)
    Fr = GF(prime_r)

    # Base field
    print('prime_q = {}'.format(prime_q))
    #params_generator.generate_libff_Fp_model_params(prime_q)
    Fq = GF(prime_q)

    # E/Fq
    curve = EllipticCurve(Fq, [0, 1])
    curve_order = curve.order()

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

    g1_t_frob = g1_t(param)
    print(f"g1_t_frob = {g1_t_frob}")

    t_minus_1_over_2 = (t-1)/2
    print('t_minus_1_over_2 = {}'.format(t_minus_1_over_2))

    # Build the twist E'/Fq2
    # Fq2 is constructed as Fq[u]/(u^2 - beta)Fq[u], where beta = -5
    fq_non_residue = Fq(-5)
    print('fq_non_residue = {}'.format(fq_non_residue))
    Fqx.<j> = PolynomialRing(Fq, 'j')
    assert(Fqx(j^2 - fq_non_residue).is_irreducible())
    Fq2.<u> = GF(prime_q^2, modulus=j^2 - fq_non_residue)

    nqr = Fq2(u)
    assert(nqr^euler_fp2 == Fq2(-1))
    nqr_to_t = pow(nqr, t, Fq2(u^2 + 5))
    print('nqr_to_t = {}'.format(nqr_to_t))

    m_twist = EllipticCurve(Fq2, [0, (1 * u)])
    d_twist = EllipticCurve(Fq2, [0, (1 / u)])
    twist = 0
    if (m_twist.order() % prime_r) == 0:
        print("Good twist is: M-type")
        twist = m_twist
    elif (d_twist.order() % prime_r) == 0:
        print("Good twist is: D-type")
        twist = d_twist
    else:
        raise BaseException("Error. None of the proposed twists has order divisible by r")

    frob_coeff_c1_1 = fq_non_residue^((0/2)) # = 1
    frob_coeff_c1_2 = fq_non_residue^(((prime_q) - 1)/2)

    twist_order = twist.order()
    assert(g1_order(curve_order) == prime_r)

    # Cofactors
    h1 = g1_h(param)
    # Check that the cofactor formula is sound
    assert(h1 == curve_order // prime_r)
    print(f"G1 order: {curve_order}")
    print('h1 = {}'.format(h1))
    assert(curve_order == prime_q + 1 - g1_t_frob)

    h2 = g2_h(param)
    # Check that the cofactor formula is sound
    assert(h2 == twist_order // prime_r)
    print(f"G2 twist order: {twist_order}")
    print('h2 = {}'.format(h2))

    # G1 fast subgroup check coefficients
    print("G1 fast subgroup check coefficients:")
    r1_0, r1_1 = g1_fast_subgroup_check_coefficients(curve_order, prime_q, prime_r)
    print(f" r1_0: {r1_0}")
    print(f" r1_1: {r1_1} ({hex(r1_1)})")

    # G2 fast cofactor multiplication coefficients
    print(f"G2 fast mul_by_cofactor coefficients:")
    h2_0, h2_1 = g2_fast_cofactor_coefficients(twist_order, prime_q, h2)
    print(f" h2_0={h2_0}")
    print(f" h2_1={h2_1}")

    # h2_0, h2_1 = fuentes_g2_fast_cofactor_coefficients(param)
    # print(f"[Fuentes] G2 fast mul_by_cofactor coefficients: h2_0={h2_0}, h2_1={h2_1}")

    # G2 fast subgroup check coefficients
    print("G2 fast subgroup check coefficients:")
    r2_0, r2_1, r2_2 = g2_fast_subgroup_check_coefficients(twist_order, prime_q, prime_r, g1_t_frob, h2)
    print(f" r2_0={r2_0}")
    print(f" r2_1={r2_1}")
    print(f" r2_2={r2_2}")

    # Build the full tower extension
    """
    # Fq6 is constructed as Fq2[v]/(v^3 - u)Fq2[v]
    fq2_non_residue = Fq2(u)
    Fq2x.<l> = PolynomialRing(Fq2, 'l')
    assert(Fq2x(l^3 - fq2_non_residue).is_irreducible())
    # Fq6.<v> = GF((prime_q^2)^3, modulus=Fq2x(l^3 - fq2_non_residue))
    fq6_non_residue = Fq6(v)
    print('non_residue = {}'.format(non_residue))

    # Fq12 is constructed as Fq6[w]/(w^2 + v)F[w]
    Fq3x.<w> = PolynomialRing(Fq3, 'w')
    assert(Fqx(w^2 - v).is_irreducible())
    Fq12.<n> = GF(((prime_q^2)^3)^2, modulus=w^2 -v
    """

if __name__ == '__main__':
    main()
