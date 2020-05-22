import generate_libff_Fp_model_params

# Curve parameter
u = 0x8508c00000000001

def t_0(x):
    return x**5 - 3*x**4 + 3*x**3 - x + 3

def y_0(x):
    return t_0(x)/3

def t_1(x):
    return -x**5 + 3*x**4 - 3*x**3 + x

def y_1(x):
    return (-t_1(x))/3

def t_i(x, i):
    if i not in {0, 1}:
        raise BaseException("i must be 0 or 1")
    if i == 0:
        return t_0(x)
    return t_1(x)

def y_i(x, i):
    if i not in {0, 1}:
        raise BaseException("i must be 0 or 1")
    if i == 0:
        return y_0(x)
    return y_1(x)

# Lifting cofactors choosen
h_t = 13
h_y = 9

# Choosen set of solutions (i \in {0,1})
i = 0

# See Table 4: https://eprint.iacr.org/2020/351.pdf
def r(x):
    return (x**6 - 2*x**5 + 2*x**3 + x + 1)/3

def t(x):
    return r(x) * h_t + t_i(x, i)

def y(x):
    return r(x) * h_y + y_i(x, i)

def q(x):
    return (t(x)**2 + 3*y(x)**2)/4

prime_r = r(u)
generate_libff_Fp_model_params(prime_r)
Fr = GF(prime_r)

prime_q = q(u)
generate_libff_Fp_model_params(prime_q)
Fq = GF(prime_q)