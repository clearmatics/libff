/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/alt_bn128/alt_bn128_g1.hpp>

namespace libff
{

static const uint8_t G1_ZERO_FLAG = 1 << 0;
static const uint8_t G1_Y_LSB_FLAG = 1 << 1;

#ifdef PROFILE_OP_COUNTS
long long alt_bn128_G1::add_cnt = 0;
long long alt_bn128_G1::dbl_cnt = 0;
#endif

std::vector<size_t> alt_bn128_G1::wnaf_window_table;
std::vector<size_t> alt_bn128_G1::fixed_base_exp_window_table;
alt_bn128_G1 alt_bn128_G1::G1_zero;
alt_bn128_G1 alt_bn128_G1::G1_one;
alt_bn128_Fq alt_bn128_G1::coeff_a;
alt_bn128_Fq alt_bn128_G1::coeff_b;
bigint<alt_bn128_G1::h_limbs> alt_bn128_G1::h;

alt_bn128_G1::alt_bn128_G1()
{
    this->X = G1_zero.X;
    this->Y = G1_zero.Y;
    this->Z = G1_zero.Z;
}

void alt_bn128_G1::print() const
{
    if (this->is_zero()) {
        printf("O\n");
    } else {
        alt_bn128_G1 copy(*this);
        copy.to_affine_coordinates();
        gmp_printf(
            "(%Nd , %Nd)\n",
            copy.X.as_bigint().data,
            alt_bn128_Fq::num_limbs,
            copy.Y.as_bigint().data,
            alt_bn128_Fq::num_limbs);
    }
}

void alt_bn128_G1::print_coordinates() const
{
    if (this->is_zero()) {
        printf("O\n");
    } else {
        gmp_printf(
            "(%Nd : %Nd : %Nd)\n",
            this->X.as_bigint().data,
            alt_bn128_Fq::num_limbs,
            this->Y.as_bigint().data,
            alt_bn128_Fq::num_limbs,
            this->Z.as_bigint().data,
            alt_bn128_Fq::num_limbs);
    }
}

void alt_bn128_G1::to_affine_coordinates()
{
    if (this->is_zero()) {
        this->X = alt_bn128_Fq::zero();
        this->Y = alt_bn128_Fq::one();
        this->Z = alt_bn128_Fq::zero();
    } else {
        const alt_bn128_Fq Z_inv = Z.inverse();
        const alt_bn128_Fq Z2_inv = Z_inv.squared();
        const alt_bn128_Fq Z3_inv = Z2_inv * Z_inv;
        this->X = this->X * Z2_inv;
        this->Y = this->Y * Z3_inv;
        this->Z = alt_bn128_Fq::one();
    }
}

void alt_bn128_G1::to_special() { this->to_affine_coordinates(); }

bool alt_bn128_G1::is_special() const
{
    return (this->is_zero() || this->Z == alt_bn128_Fq::one());
}

bool alt_bn128_G1::is_zero() const { return (this->Z.is_zero()); }

bool alt_bn128_G1::operator==(const alt_bn128_G1 &other) const
{
    if (this->is_zero()) {
        return other.is_zero();
    }

    if (other.is_zero()) {
        return false;
    }

    // now neither is O

    // using Jacobian coordinates so:
    //   (X1:Y1:Z1) = (X2:Y2:Z2)
    //   iff
    //   X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    //   iff
    //   X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    const alt_bn128_Fq Z1_squared = (this->Z).squared();
    const alt_bn128_Fq Z2_squared = (other.Z).squared();

    if ((this->X * Z2_squared) != (other.X * Z1_squared)) {
        return false;
    }

    const alt_bn128_Fq Z1_cubed = (this->Z) * Z1_squared;
    const alt_bn128_Fq Z2_cubed = (other.Z) * Z2_squared;

    if ((this->Y * Z2_cubed) != (other.Y * Z1_cubed)) {
        return false;
    }

    return true;
}

bool alt_bn128_G1::operator!=(const alt_bn128_G1 &other) const
{
    return !(operator==(other));
}

alt_bn128_G1 alt_bn128_G1::operator+(const alt_bn128_G1 &other) const
{
    return add(other);
}

alt_bn128_G1 alt_bn128_G1::operator-() const
{
    return alt_bn128_G1(this->X, -(this->Y), this->Z);
}

alt_bn128_G1 alt_bn128_G1::operator-(const alt_bn128_G1 &other) const
{
    return (*this) + (-other);
}

alt_bn128_G1 alt_bn128_G1::add(const alt_bn128_G1 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero()) {
        return other;
    }

    if (other.is_zero()) {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // handle double case
    if (this->operator==(other)) {
        return this->dbl();
    }

#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl

    // Z1Z1 = Z1^2
    const alt_bn128_Fq Z1Z1 = (this->Z).squared();
    // Z2Z2 = Z2^2
    const alt_bn128_Fq Z2Z2 = (other.Z).squared();
    // U1 = X1 * Z2Z2
    const alt_bn128_Fq U1 = (this->X) * Z2Z2;
    // U2 = X2 * Z1Z1
    const alt_bn128_Fq U2 = (other.X) * Z1Z1;
    // S1 = Y1 * Z2 * Z2Z2
    const alt_bn128_Fq S1 = (this->Y) * (other.Z) * Z2Z2;
    // S2 = Y2 * Z1 * Z1Z1
    const alt_bn128_Fq S2 = (other.Y) * (this->Z) * Z1Z1;
    // H = U2-U1
    const alt_bn128_Fq H = U2 - U1;
    const alt_bn128_Fq S2_minus_S1 = S2 - S1;
    // I = (2 * H)^2
    const alt_bn128_Fq I = (H + H).squared();
    // J = H * I
    const alt_bn128_Fq J = H * I;
    // r = 2 * (S2-S1)
    const alt_bn128_Fq r = S2_minus_S1 + S2_minus_S1;
    // V = U1 * I
    const alt_bn128_Fq V = U1 * I;
    // X3 = r^2 - J - 2 * V
    const alt_bn128_Fq X3 = r.squared() - J - (V + V);
    const alt_bn128_Fq S1_J = S1 * J;
    // Y3 = r * (V-X3)-2 S1 J
    const alt_bn128_Fq Y3 = r * (V - X3) - (S1_J + S1_J);
    // Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2) * H
    const alt_bn128_Fq Z3 = ((this->Z + other.Z).squared() - Z1Z1 - Z2Z2) * H;

    return alt_bn128_G1(X3, Y3, Z3);
}

alt_bn128_G1 alt_bn128_G1::mixed_add(const alt_bn128_G1 &other) const
{
#ifdef DEBUG
    assert(other.is_special());
#endif

    // handle special cases having to do with O
    if (this->is_zero()) {
        return other;
    }

    if (other.is_zero()) {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // check for doubling case

    // using Jacobian coordinates so:
    //   (X1:Y1:Z1) = (X2:Y2:Z2)
    //   iff
    //   X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    //   iff
    //   X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    // we know that Z2 = 1

    const alt_bn128_Fq Z1Z1 = (this->Z).squared();

    const alt_bn128_Fq &U1 = this->X;
    const alt_bn128_Fq U2 = other.X * Z1Z1;

    const alt_bn128_Fq Z1_cubed = (this->Z) * Z1Z1;

    // S1 = Y1 * Z2 * Z2Z2
    const alt_bn128_Fq &S1 = (this->Y);
    // S2 = Y2 * Z1 * Z1Z1
    const alt_bn128_Fq S2 = (other.Y) * Z1_cubed;

    if (U1 == U2 && S1 == S2) {
        // dbl case; nothing of above can be reused
        return this->dbl();
    }

#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif

    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
    // H = U2-X1
    const alt_bn128_Fq H = U2 - (this->X);
    // HH = H&2
    const alt_bn128_Fq HH = H.squared();
    // I = 4*HH
    alt_bn128_Fq I = HH + HH;
    I = I + I;
    // J = H*I
    const alt_bn128_Fq J = H * I;
    // r = 2*(S2-Y1)
    alt_bn128_Fq r = S2 - (this->Y);
    r = r + r;
    // V = X1*I
    const alt_bn128_Fq V = (this->X) * I;
    // X3 = r^2-J-2*V
    const alt_bn128_Fq X3 = r.squared() - J - V - V;
    // Y3 = r*(V-X3)-2*Y1*J
    alt_bn128_Fq Y3 = (this->Y) * J;
    Y3 = r * (V - X3) - Y3 - Y3;
    // Z3 = (Z1+H)^2-Z1Z1-HH
    const alt_bn128_Fq Z3 = ((this->Z) + H).squared() - Z1Z1 - HH;

    return alt_bn128_G1(X3, Y3, Z3);
}

alt_bn128_G1 alt_bn128_G1::dbl() const
{
#ifdef PROFILE_OP_COUNTS
    this->dbl_cnt++;
#endif
    // handle point at infinity
    if (this->is_zero()) {
        return (*this);
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l

    // A = X1^2
    const alt_bn128_Fq A = (this->X).squared();
    // B = Y1^2
    const alt_bn128_Fq B = (this->Y).squared();
    // C = B^2
    const alt_bn128_Fq C = B.squared();
    alt_bn128_Fq D = (this->X + B).squared() - A - C;
    // D = 2 * ((X1 + B)^2 - A - C)
    D = D + D;
    // E = 3 * A
    const alt_bn128_Fq E = A + A + A;
    // F = E^2
    const alt_bn128_Fq F = E.squared();
    // X3 = F - 2 D
    const alt_bn128_Fq X3 = F - (D + D);
    alt_bn128_Fq eightC = C + C;
    eightC = eightC + eightC;
    eightC = eightC + eightC;
    // Y3 = E * (D - X3) - 8 * C
    const alt_bn128_Fq Y3 = E * (D - X3) - eightC;
    const alt_bn128_Fq Y1Z1 = (this->Y) * (this->Z);
    // Z3 = 2 * Y1 * Z1
    const alt_bn128_Fq Z3 = Y1Z1 + Y1Z1;

    return alt_bn128_G1(X3, Y3, Z3);
}

alt_bn128_G1 alt_bn128_G1::mul_by_cofactor() const
{
    // Cofactor = 1
    return *this;
}

bool alt_bn128_G1::is_well_formed() const
{
    if (this->is_zero()) {
        return true;
    } else {
        // y^2 = x^3 + b
        //
        // We are using Jacobian coordinates, so equation we need to check is
        // actually
        //
        // (y/z^3)^2 = (x/z^2)^3 + b
        // y^2 / z^6 = x^3 / z^6 + b
        // y^2 = x^3 + b z^6
        const alt_bn128_Fq X2 = this->X.squared();
        const alt_bn128_Fq Y2 = this->Y.squared();
        const alt_bn128_Fq Z2 = this->Z.squared();

        const alt_bn128_Fq X3 = this->X * X2;
        const alt_bn128_Fq Z3 = this->Z * Z2;
        const alt_bn128_Fq Z6 = Z3.squared();

        return (Y2 == X3 + alt_bn128_coeff_b * Z6);
    }
}

bool alt_bn128_G1::is_in_safe_subgroup() const
{
    // G1 is the entire group E(Fq), and so there is nothing to check here.
    return true;
}

const alt_bn128_G1 &alt_bn128_G1::zero() { return G1_zero; }

const alt_bn128_G1 &alt_bn128_G1::one() { return G1_one; }

alt_bn128_G1 alt_bn128_G1::random_element()
{
    return (scalar_field::random_element().as_bigint()) * G1_one;
}

void alt_bn128_G1::write_uncompressed(std::ostream &out) const
{
    // <is_zero> | <x-coord> | <y-coord>
    alt_bn128_G1 copy((*this));
    copy.to_affine_coordinates();
    const char is_zero = '0' + (copy.is_zero() ? G1_ZERO_FLAG : 0);
    out.write(&is_zero, 1) << copy.X << copy.Y;
}

void alt_bn128_G1::write_compressed(std::ostream &out) const
{
    // <flags> | <x-coord>
    alt_bn128_G1 copy((*this));
    copy.to_affine_coordinates();
    const uint8_t flags =
        (copy.is_zero() ? G1_ZERO_FLAG : 0) |
        ((copy.Y.as_bigint().data[0] & 1) ? G1_Y_LSB_FLAG : 0);
    const char flags_char = '0' + flags;
    out.write(&flags_char, 1) << copy.X;
}

void alt_bn128_G1::read_uncompressed(std::istream &in, alt_bn128_G1 &out)
{
    // <is_zero> | <x-coord> | <y-coord>
    char is_zero_char;
    in.read(&is_zero_char, 1) >> out.X >> out.Y;
    const uint8_t is_zero = (is_zero_char - '0') & G1_ZERO_FLAG;

    if (!is_zero) {
        out.Z = alt_bn128_Fq::one();
    } else {
        out = alt_bn128_G1::zero();
    }
}

void alt_bn128_G1::read_compressed(std::istream &in, alt_bn128_G1 &out)
{
    // <flags> | <x-coord>
    char flags_char;
    in.read(&flags_char, 1) >> out.X;
    ;
    consume_OUTPUT_SEPARATOR(in);
    const uint8_t flags = flags_char - '0';

    // y = +/- sqrt(x^3 + b)
    if (0 == (flags & G1_ZERO_FLAG)) {
        // flags >> 1
        const uint8_t Y_lsb = 0 != (flags & G1_Y_LSB_FLAG);
        const alt_bn128_Fq tX2 = out.X.squared();
        const alt_bn128_Fq tY2 = tX2 * out.X + alt_bn128_coeff_b;
        out.Y = tY2.sqrt();

        if ((uint8_t)(out.Y.as_bigint().data[0] & 1) != Y_lsb) {
            out.Y = -out.Y;
        }

        out.Z = alt_bn128_Fq::one();
    } else {
        out = alt_bn128_G1::zero();
    }
}

std::ostream &operator<<(std::ostream &out, const alt_bn128_G1 &g)
{
#ifdef NO_PT_COMPRESSION
    g.write_uncompressed(out);
#else
    g.write_compressed(out);
#endif
    return out;
}

std::istream &operator>>(std::istream &in, alt_bn128_G1 &g)
{
#ifdef NO_PT_COMPRESSION
    alt_bn128_G1::read_uncompressed(in, g);
#else
    alt_bn128_G1::read_compressed(in, g);
#endif
    return in;
}

void alt_bn128_G1::batch_to_special_all_non_zeros(
    std::vector<alt_bn128_G1> &vec)
{
    std::vector<alt_bn128_Fq> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el : vec) {
        Z_vec.emplace_back(el.Z);
    }
    batch_invert<alt_bn128_Fq>(Z_vec);

    const alt_bn128_Fq one = alt_bn128_Fq::one();

    for (size_t i = 0; i < vec.size(); ++i) {
        const alt_bn128_Fq Z2 = Z_vec[i].squared();
        const alt_bn128_Fq Z3 = Z_vec[i] * Z2;

        vec[i].X = vec[i].X * Z2;
        vec[i].Y = vec[i].Y * Z3;
        vec[i].Z = one;
    }
}

} // namespace libff
