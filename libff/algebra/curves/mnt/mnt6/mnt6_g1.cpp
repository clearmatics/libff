/** @file
 *****************************************************************************

 Implementation of interfaces for the MNT6 G1 group.

 See mnt6_g1.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/mnt/mnt6/mnt6_g1.hpp>

namespace libff
{

#ifdef PROFILE_OP_COUNTS
long long mnt6_G1::add_cnt = 0;
long long mnt6_G1::dbl_cnt = 0;
#endif

std::vector<size_t> mnt6_G1::wnaf_window_table;
std::vector<size_t> mnt6_G1::fixed_base_exp_window_table;
mnt6_G1 mnt6_G1::G1_zero;
mnt6_G1 mnt6_G1::G1_one;
mnt6_Fq mnt6_G1::coeff_a;
mnt6_Fq mnt6_G1::coeff_b;
bigint<mnt6_G1::h_limbs> mnt6_G1::h;

mnt6_G1::mnt6_G1()
{
    this->X = G1_zero.X;
    this->Y = G1_zero.Y;
    this->Z = G1_zero.Z;
}

void mnt6_G1::print() const
{
    if (this->is_zero()) {
        printf("O\n");
    } else {
        mnt6_G1 copy(*this);
        copy.to_affine_coordinates();
        gmp_printf(
            "(%Nd , %Nd)\n",
            copy.X.as_bigint().data,
            mnt6_Fq::num_limbs,
            copy.Y.as_bigint().data,
            mnt6_Fq::num_limbs);
    }
}

void mnt6_G1::print_coordinates() const
{
    if (this->is_zero()) {
        printf("O\n");
    } else {
        gmp_printf(
            "(%Nd : %Nd : %Nd)\n",
            this->X.as_bigint().data,
            mnt6_Fq::num_limbs,
            this->Y.as_bigint().data,
            mnt6_Fq::num_limbs,
            this->Z.as_bigint().data,
            mnt6_Fq::num_limbs);
    }
}

void mnt6_G1::to_affine_coordinates()
{
    if (this->is_zero()) {
        this->X = mnt6_Fq::zero();
        this->Y = mnt6_Fq::one();
        this->Z = mnt6_Fq::zero();
    } else {
        const mnt6_Fq Z_inv = Z.inverse();
        this->X = this->X * Z_inv;
        this->Y = this->Y * Z_inv;
        this->Z = mnt6_Fq::one();
    }
}

void mnt6_G1::to_special() { this->to_affine_coordinates(); }

bool mnt6_G1::is_special() const
{
    return (this->is_zero() || this->Z == mnt6_Fq::one());
}

bool mnt6_G1::is_zero() const
{
    return (this->X.is_zero() && this->Z.is_zero());
}

bool mnt6_G1::operator==(const mnt6_G1 &other) const
{
    if (this->is_zero()) {
        return other.is_zero();
    }

    if (other.is_zero()) {
        return false;
    }

    /* now neither is O */

    // X1/Z1 = X2/Z2 <=> X1*Z2 = X2*Z1
    if ((this->X * other.Z) != (other.X * this->Z)) {
        return false;
    }

    // Y1/Z1 = Y2/Z2 <=> Y1*Z2 = Y2*Z1
    if ((this->Y * other.Z) != (other.Y * this->Z)) {
        return false;
    }

    return true;
}

bool mnt6_G1::operator!=(const mnt6_G1 &other) const
{
    return !(operator==(other));
}

mnt6_G1 mnt6_G1::operator+(const mnt6_G1 &other) const
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

    // handle double case, and then all the rest

    // The code below is equivalent to (but faster than) the snippet below:
    //
    // if (this->operator==(other))
    // {
    //   return this->dbl();
    // }
    // else
    // {
    //   return this->add(other);
    // }

    // X1Z2 = X1*Z2
    const mnt6_Fq X1Z2 = (this->X) * (other.Z);
    // X2Z1 = X2*Z1
    const mnt6_Fq X2Z1 = (this->Z) * (other.X);

    // (used both in add and double checks)

    // Y1Z2 = Y1*Z2
    const mnt6_Fq Y1Z2 = (this->Y) * (other.Z);
    // Y2Z1 = Y2*Z1
    const mnt6_Fq Y2Z1 = (this->Z) * (other.Y);

    if (X1Z2 == X2Z1 && Y1Z2 == Y2Z1) {
        // perform dbl case
        // XX  = X1^2
        const mnt6_Fq XX = (this->X).squared();
        // ZZ  = Z1^2
        const mnt6_Fq ZZ = (this->Z).squared();
        // w   = a*ZZ + 3*XX
        const mnt6_Fq w = mnt6_G1::coeff_a * ZZ + (XX + XX + XX);
        const mnt6_Fq Y1Z1 = (this->Y) * (this->Z);
        // s   = 2*Y1*Z1
        const mnt6_Fq s = Y1Z1 + Y1Z1;
        // ss  = s^2
        const mnt6_Fq ss = s.squared();
        // sss = s*ss
        const mnt6_Fq sss = s * ss;
        // R   = Y1*s
        const mnt6_Fq R = (this->Y) * s;
        // RR  = R^2
        const mnt6_Fq RR = R.squared();
        // B   = (X1+R)^2 - XX - RR
        const mnt6_Fq B = ((this->X) + R).squared() - XX - RR;
        // h   = w^2 - 2*B
        const mnt6_Fq h = w.squared() - (B + B);
        // X3  = h*s
        const mnt6_Fq X3 = h * s;
        // Y3  = w*(B-h) - 2*RR
        const mnt6_Fq Y3 = w * (B - h) - (RR + RR);
        // Z3  = sss
        const mnt6_Fq Z3 = sss;

        return mnt6_G1(X3, Y3, Z3);
    }

    // if we have arrived here we are in the add case
    // Z1Z2 = Z1*Z2
    const mnt6_Fq Z1Z2 = (this->Z) * (other.Z);
    // u    = Y2*Z1-Y1Z2
    const mnt6_Fq u = Y2Z1 - Y1Z2;
    // uu   = u^2
    const mnt6_Fq uu = u.squared();
    // v    = X2*Z1-X1Z2
    const mnt6_Fq v = X2Z1 - X1Z2;
    // vv   = v^2
    const mnt6_Fq vv = v.squared();
    // vvv  = v*vv
    const mnt6_Fq vvv = v * vv;
    // R    = vv*X1Z2
    const mnt6_Fq R = vv * X1Z2;
    // A    = uu*Z1Z2 - vvv - 2*R
    const mnt6_Fq A = uu * Z1Z2 - (vvv + R + R);
    // X3   = v*A
    const mnt6_Fq X3 = v * A;
    // Y3   = u*(R-A) - vvv*Y1Z2
    const mnt6_Fq Y3 = u * (R - A) - vvv * Y1Z2;
    // Z3   = vvv*Z1Z2
    const mnt6_Fq Z3 = vvv * Z1Z2;

    return mnt6_G1(X3, Y3, Z3);
}

mnt6_G1 mnt6_G1::operator-() const
{
    return mnt6_G1(this->X, -(this->Y), this->Z);
}

mnt6_G1 mnt6_G1::operator-(const mnt6_G1 &other) const
{
    return (*this) + (-other);
}

mnt6_G1 mnt6_G1::add(const mnt6_G1 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero()) {
        return other;
    }

    if (other.is_zero()) {
        return (*this);
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
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-add-1998-cmo-2

    // Y1Z2 = Y1*Z2
    const mnt6_Fq Y1Z2 = (this->Y) * (other.Z);
    // X1Z2 = X1*Z2
    const mnt6_Fq X1Z2 = (this->X) * (other.Z);
    // Z1Z2 = Z1*Z2
    const mnt6_Fq Z1Z2 = (this->Z) * (other.Z);
    // u    = Y2*Z1-Y1Z2
    const mnt6_Fq u = (other.Y) * (this->Z) - Y1Z2;
    // uu   = u^2
    const mnt6_Fq uu = u.squared();
    // v    = X2*Z1-X1Z2
    const mnt6_Fq v = (other.X) * (this->Z) - X1Z2;
    // vv   = v^2
    const mnt6_Fq vv = v.squared();
    // vvv  = v*vv
    const mnt6_Fq vvv = v * vv;
    // R    = vv*X1Z2
    const mnt6_Fq R = vv * X1Z2;
    // A    = uu*Z1Z2 - vvv - 2*R
    const mnt6_Fq A = uu * Z1Z2 - (vvv + R + R);
    // X3   = v*A
    const mnt6_Fq X3 = v * A;
    // Y3   = u*(R-A) - vvv*Y1Z2
    const mnt6_Fq Y3 = u * (R - A) - vvv * Y1Z2;
    // Z3   = vvv*Z1Z2
    const mnt6_Fq Z3 = vvv * Z1Z2;

    return mnt6_G1(X3, Y3, Z3);
}

mnt6_G1 mnt6_G1::mixed_add(const mnt6_G1 &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-add-1998-cmo-2
    // assert(other.Z == mnt6_Fq::one());

    if (this->is_zero()) {
        return other;
    }

    if (other.is_zero()) {
        return (*this);
    }

#ifdef DEBUG
    assert(other.is_special());
#endif

    // X1Z2 = X1*Z2 (but other is special and not zero)
    const mnt6_Fq &X1Z2 = (this->X);
    // X2Z1 = X2*Z1
    const mnt6_Fq X2Z1 = (this->Z) * (other.X);

    // (used both in add and double checks)

    // Y1Z2 = Y1*Z2 (but other is special and not zero)
    const mnt6_Fq &Y1Z2 = (this->Y);
    // Y2Z1 = Y2*Z1
    const mnt6_Fq Y2Z1 = (this->Z) * (other.Y);

    if (X1Z2 == X2Z1 && Y1Z2 == Y2Z1) {
        return this->dbl();
    }

    // u = Y2*Z1-Y1
    mnt6_Fq u = Y2Z1 - this->Y;
    // uu = u2
    mnt6_Fq uu = u.squared();
    // v = X2*Z1-X1
    mnt6_Fq v = X2Z1 - this->X;
    // vv = v2
    mnt6_Fq vv = v.squared();
    // vvv = v*vv
    mnt6_Fq vvv = v * vv;
    // R = vv*X1
    mnt6_Fq R = vv * this->X;
    // A = uu*Z1-vvv-2*R
    mnt6_Fq A = uu * this->Z - vvv - R - R;
    // X3 = v*A
    mnt6_Fq X3 = v * A;
    // Y3 = u*(R-A)-vvv*Y1
    mnt6_Fq Y3 = u * (R - A) - vvv * this->Y;
    // Z3 = vvv*Z1
    mnt6_Fq Z3 = vvv * this->Z;

    return mnt6_G1(X3, Y3, Z3);
}

mnt6_G1 mnt6_G1::dbl() const
{
#ifdef PROFILE_OP_COUNTS
    this->dbl_cnt++;
#endif
    if (this->is_zero()) {
        return (*this);
    } else {
        // NOTE: does not handle O and pts of order 2,4
        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#doubling-dbl-2007-bl

        // XX  = X1^2
        const mnt6_Fq XX = (this->X).squared();
        // ZZ  = Z1^2
        const mnt6_Fq ZZ = (this->Z).squared();
        // w   = a*ZZ + 3*XX
        const mnt6_Fq w = mnt6_G1::coeff_a * ZZ + (XX + XX + XX);
        const mnt6_Fq Y1Z1 = (this->Y) * (this->Z);
        // s   = 2*Y1*Z1
        const mnt6_Fq s = Y1Z1 + Y1Z1;
        // ss  = s^2
        const mnt6_Fq ss = s.squared();
        // sss = s*ss
        const mnt6_Fq sss = s * ss;
        // R   = Y1*s
        const mnt6_Fq R = (this->Y) * s;
        // RR  = R^2
        const mnt6_Fq RR = R.squared();
        // B   = (X1+R)^2 - XX - RR
        const mnt6_Fq B = ((this->X) + R).squared() - XX - RR;
        // h   = w^2 - 2*B
        const mnt6_Fq h = w.squared() - (B + B);
        // X3  = h*s
        const mnt6_Fq X3 = h * s;
        // Y3  = w*(B-h) - 2*RR
        const mnt6_Fq Y3 = w * (B - h) - (RR + RR);
        // Z3  = sss
        const mnt6_Fq Z3 = sss;

        return mnt6_G1(X3, Y3, Z3);
    }
}

mnt6_G1 mnt6_G1::mul_by_cofactor() const
{
    // Cofactor = 1
    return (*this);
}

bool mnt6_G1::is_well_formed() const
{
    if (this->is_zero()) {
        return true;
    } else {
        // y^2 = x^3 + ax + b
        //
        // We are using projective, so equation we need to check is actually
        //
        // (y/z)^2 = (x/z)^3 + a (x/z) + b
        // z y^2 = x^3  + a z^2 x + b z^3
        //
        // z (y^2 - b z^2) = x ( x^2 + a z^2)
        const mnt6_Fq X2 = this->X.squared();
        const mnt6_Fq Y2 = this->Y.squared();
        const mnt6_Fq Z2 = this->Z.squared();

        return (
            this->Z * (Y2 - mnt6_G1::coeff_b * Z2) ==
            this->X * (X2 + mnt6_G1::coeff_a * Z2));
    }
}

bool mnt6_G1::is_in_safe_subgroup() const { return true; }

const mnt6_G1 &mnt6_G1::zero() { return G1_zero; }

const mnt6_G1 &mnt6_G1::one() { return G1_one; }

mnt6_G1 mnt6_G1::random_element()
{
    return (scalar_field::random_element().as_bigint()) * G1_one;
}

void mnt6_G1::write_uncompressed(std::ostream &out) const
{
    mnt6_G1 copy(*this);
    copy.to_affine_coordinates();

    out << (copy.is_zero() ? 1 : 0) << OUTPUT_SEPARATOR;
    out << copy.X << OUTPUT_SEPARATOR << copy.Y;
}

void mnt6_G1::write_compressed(std::ostream &out) const
{
    mnt6_G1 copy(*this);
    copy.to_affine_coordinates();

    out << (copy.is_zero() ? 1 : 0) << OUTPUT_SEPARATOR;
    /* storing LSB of Y */
    out << copy.X << OUTPUT_SEPARATOR << (copy.Y.as_bigint().data[0] & 1);
}

void mnt6_G1::read_uncompressed(std::istream &in, mnt6_G1 &g)
{
    char is_zero;
    mnt6_Fq tX, tY;
    in >> is_zero >> tX >> tY;
    is_zero -= '0';

    // using projective coordinates
    if (!is_zero) {
        g.X = tX;
        g.Y = tY;
        g.Z = mnt6_Fq::one();
    } else {
        g = mnt6_G1::zero();
    }
}

void mnt6_G1::read_compressed(std::istream &in, mnt6_G1 &g)
{
    char is_zero;
    mnt6_Fq tX, tY;
    // this reads is_zero;
    in.read((char *)&is_zero, 1);
    is_zero -= '0';
    consume_OUTPUT_SEPARATOR(in);

    unsigned char Y_lsb;
    in >> tX;
    consume_OUTPUT_SEPARATOR(in);
    in.read((char *)&Y_lsb, 1);
    Y_lsb -= '0';

    // y = +/- sqrt(x^3 + a*x + b)
    if (!is_zero) {
        mnt6_Fq tX2 = tX.squared();
        mnt6_Fq tY2 = (tX2 + mnt6_G1::coeff_a) * tX + mnt6_G1::coeff_b;
        tY = tY2.sqrt();

        if ((tY.as_bigint().data[0] & 1) != Y_lsb) {
            tY = -tY;
        }
    }

    // using projective coordinates
    if (!is_zero) {
        g.X = tX;
        g.Y = tY;
        g.Z = mnt6_Fq::one();
    } else {
        g = mnt6_G1::zero();
    }
}

void mnt6_G1::batch_to_special_all_non_zeros(std::vector<mnt6_G1> &vec)
{
    std::vector<mnt6_Fq> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el : vec) {
        Z_vec.emplace_back(el.Z);
    }
    batch_invert<mnt6_Fq>(Z_vec);

    const mnt6_Fq one = mnt6_Fq::one();

    for (size_t i = 0; i < vec.size(); ++i) {
        vec[i] = mnt6_G1(vec[i].X * Z_vec[i], vec[i].Y * Z_vec[i], one);
    }
}

std::ostream &operator<<(std::ostream &out, const mnt6_G1 &g)
{
#ifdef NO_PT_COMPRESSION
    g.write_uncompressed(out);
#else
    g.write_compressed(out);
#endif
    return out;
}

std::istream &operator>>(std::istream &in, mnt6_G1 &g)
{
#ifdef NO_PT_COMPRESSION
    mnt6_G1::read_uncompressed(in, g);
#else
    mnt6_G1::read_compressed(in, g);
#endif
    return in;
}

} // namespace libff
