#include <libff/algebra/curves/bls12_377/bls12_377_g1.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_g2.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_init.hpp>

// Note: These parameters match the RUST implementation of the BLS12-377 curve:
// https://github.com/scipr-lab/zexe/tree/6bfe574f7adea14b97ff554bbb594988635b1908/algebra/src/bls12_377

namespace libff {

bigint<bls12_377_r_limbs> bls12_377_modulus_r;
// bls12_377_modulus_q is a macro referring to bw6_761_modulus_r. See
// bls12_377_init.hpp.
// bigint<bls12_377_q_limbs> bls12_377_modulus_q;

bls12_377_Fq bls12_377_coeff_b;
bigint<bls12_377_r_limbs> bls12_377_trace_of_frobenius;
bls12_377_Fq2 bls12_377_twist;
bls12_377_Fq2 bls12_377_twist_coeff_b;
bls12_377_Fq bls12_377_twist_mul_by_b_c0;
bls12_377_Fq bls12_377_twist_mul_by_b_c1;
bls12_377_Fq2 bls12_377_twist_mul_by_q_X;
bls12_377_Fq2 bls12_377_twist_mul_by_q_Y;

// See bls12_377_G1::is_in_safe_subgroup
bls12_377_Fq bls12_377_g1_endomorphism_beta;
bigint<bls12_377_r_limbs> bls12_377_g1_safe_subgroup_check_c1;
bigint<bls12_377_r_limbs> bls12_377_g1_proof_of_safe_subgroup_w;
bls12_377_Fq bls12_377_g1_proof_of_safe_subgroup_non_member_x;
bls12_377_Fq bls12_377_g1_proof_of_safe_subgroup_non_member_y;

// Coefficients for G2 untwist-frobenius-twist
bls12_377_Fq12 bls12_377_g2_untwist_frobenius_twist_w;
bls12_377_Fq12 bls12_377_g2_untwist_frobenius_twist_v;
bls12_377_Fq12 bls12_377_g2_untwist_frobenius_twist_w_3;
bls12_377_Fq12 bls12_377_g2_untwist_frobenius_twist_v_inverse;
bls12_377_Fq12 bls12_377_g2_untwist_frobenius_twist_w_3_inverse;

// Coefficients used in bls12_377_G2::mul_by_cofactor
bigint<bls12_377_r_limbs> bls12_377_g2_mul_by_cofactor_h2_0;
bigint<bls12_377_r_limbs> bls12_377_g2_mul_by_cofactor_h2_1;

bigint<bls12_377_q_limbs> bls12_377_ate_loop_count;
bool bls12_377_ate_is_loop_count_neg;
// k (embedding degree) = 12
bigint<12*bls12_377_q_limbs> bls12_377_final_exponent;
bigint<bls12_377_q_limbs> bls12_377_final_exponent_z;
bool bls12_377_final_exponent_is_z_neg;

void init_bls12_377_params()
{
    typedef bigint<bls12_377_r_limbs> bigint_r;
    typedef bigint<bls12_377_q_limbs> bigint_q;

    // Montgomery assumes this
    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4);

    // Parameters for scalar field Fr
    // r = 0x12ab655e9a2ca55660b44d1e5c37b00159aa76fed00000010a11800000000001
    bls12_377_modulus_r = bigint_r("8444461749428370424248824938781546531375899335154063827935233455917409239041");
    assert(bls12_377_Fr::modulus_is_valid());
    // 64-bit architecture
    if (sizeof(mp_limb_t) == 8)
    {
        bls12_377_Fr::Rsquared = bigint_r("508595941311779472113692600146818027278633330499214071737745792929336755579");
        bls12_377_Fr::Rcubed = bigint_r("2717187485423313556320207871216538426353201097398909639086937135091399607628");
        bls12_377_Fr::inv = 0xa117fffffffffff;
    }
    // 32-bit architecture
    if (sizeof(mp_limb_t) == 4)
    {
        bls12_377_Fr::Rsquared = bigint_r("508595941311779472113692600146818027278633330499214071737745792929336755579");
        bls12_377_Fr::Rcubed = bigint_r("2717187485423313556320207871216538426353201097398909639086937135091399607628");
        bls12_377_Fr::inv = 0xffffffff;
    }
    bls12_377_Fr::num_bits = 253;
    bls12_377_Fr::euler = bigint_r("4222230874714185212124412469390773265687949667577031913967616727958704619520");
    bls12_377_Fr::s = 47;
    bls12_377_Fr::t = bigint_r("60001509534603559531609739528203892656505753216962260608619555");
    bls12_377_Fr::t_minus_1_over_2 = bigint_r("30000754767301779765804869764101946328252876608481130304309777");
    bls12_377_Fr::multiplicative_generator = bls12_377_Fr("22");
    bls12_377_Fr::root_of_unity = bls12_377_Fr("8065159656716812877374967518403273466521432693661810619979959746626482506078");
    bls12_377_Fr::nqr = bls12_377_Fr("11");
    bls12_377_Fr::nqr_to_t = bls12_377_Fr("6924886788847882060123066508223519077232160750698452411071850219367055984476");
    bls12_377_Fr::static_init();

    // Parameters for base field Fq
    // q = 0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c00000000001
    // sage: mod(0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c00000000001, 6) # = 1
    bls12_377_modulus_q = bigint_q("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177");
    assert(bls12_377_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        bls12_377_Fq::Rsquared = bigint_q("66127428376872697816332570116866232405230528984664918319606315420233909940404532140033099444330447428417853902114");
        bls12_377_Fq::Rcubed = bigint_q("157734475176213061358192738313701451942220138363611391489992831740412033225490229541667992423878570205050777755168");
        bls12_377_Fq::inv = 0x8508bfffffffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        bls12_377_Fq::Rsquared = bigint_q("66127428376872697816332570116866232405230528984664918319606315420233909940404532140033099444330447428417853902114");
        bls12_377_Fq::Rcubed = bigint_q("157734475176213061358192738313701451942220138363611391489992831740412033225490229541667992423878570205050777755168");
        bls12_377_Fq::inv = 0xffffffff;
    }

    bls12_377_Fq::num_bits = 377;
    bls12_377_Fq::euler = bigint_q("129332213006484547005326366847446766768196756377457330269942131333360234174170411387484444069786680062220160729088");
    bls12_377_Fq::s = 46;
    bls12_377_Fq::t = bigint_q("3675842578061421676390135839012792950148785745837396071634149488243117337281387659330802195819009059");
    bls12_377_Fq::t_minus_1_over_2 = bigint_q("1837921289030710838195067919506396475074392872918698035817074744121558668640693829665401097909504529");
    bls12_377_Fq::multiplicative_generator = bls12_377_Fq("15");
    bls12_377_Fq::root_of_unity = bls12_377_Fq("32863578547254505029601261939868325669770508939375122462904745766352256812585773382134936404344547323199885654433");
    // We need to find a qnr (small preferably) in order to compute square roots in the field
    bls12_377_Fq::nqr = bls12_377_Fq("5");
    bls12_377_Fq::nqr_to_t = bls12_377_Fq("33774956008227656219775876656288133547078610493828613777258829345740556592044969439504850374928261397247202212840");
    bls12_377_Fq::static_init();

    // Parameters for twist field Fq2
    bls12_377_Fq2::euler = bigint<2*bls12_377_q_limbs>("33453642642309381258089625946249069288005760010886479253070957453297957116339370141113413635838485065209570299254148838549585056123015878375022724998041828785227090063466658233059433323033772513321990316560167027213559780081664");
    bls12_377_Fq2::s = 47;
    bls12_377_Fq2::t = bigint<2*bls12_377_q_limbs>("475404855284145089315325463221726483993816145966867441829193658311651761271425728823393990805904040047516478740222806302278755994777496288961383541476974255391881599499962735436887347234371823579436839914935817251");
    bls12_377_Fq2::t_minus_1_over_2 = bigint<2*bls12_377_q_limbs>("237702427642072544657662731610863241996908072983433720914596829155825880635712864411696995402952020023758239370111403151139377997388748144480691770738487127695940799749981367718443673617185911789718419957467908625");
    // https://github.com/scipr-lab/zexe/blob/6bfe574f7adea14b97ff554bbb594988635b1908/algebra/src/bls12_377/fields/fq2.rs#L11
    // Additive inverse of 5 in GF(q)
    // sage: GF(q)(-5)
    // Fp2 = Fp[X] / (X^2 - (-5)))
    bls12_377_Fq2::non_residue = bls12_377_Fq
    ("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458172");
    bls12_377_Fq2::nqr = bls12_377_Fq2(bls12_377_Fq("0"),bls12_377_Fq("1"));
    bls12_377_Fq2::nqr_to_t = bls12_377_Fq2(bls12_377_Fq("0"),bls12_377_Fq("257286236321774568987262729980034669694531728092793737444525294935421142460394028155736019924956637466133519652786"));
    bls12_377_Fq2::Frobenius_coeffs_c1[0] = bls12_377_Fq("1");
    bls12_377_Fq2::Frobenius_coeffs_c1[1] = bls12_377_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458176");

    // Parameters for Fq6 = (Fq2)^3
    bls12_377_Fq6::non_residue = bls12_377_Fq2(bls12_377_Fq("0"),bls12_377_Fq("1"));
    bls12_377_Fq6::Frobenius_coeffs_c1[0] = bls12_377_Fq2(bls12_377_Fq("1"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c1[1] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410946"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c1[2] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c1[3] = bls12_377_Fq2(bls12_377_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458176"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c1[4] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c1[5] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047232"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[0] = bls12_377_Fq2(bls12_377_Fq("1"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[1] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[2] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[3] = bls12_377_Fq2(bls12_377_Fq("1"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[4] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[5] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),bls12_377_Fq("0"));

    // Parameters for Fq12 = ((Fq2)^3)^2
    bls12_377_Fq12::non_residue = bls12_377_Fq2(bls12_377_Fq("0"),bls12_377_Fq("1"));
    bls12_377_Fq12::Frobenius_coeffs_c1[0] = bls12_377_Fq2(bls12_377_Fq("1"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[1] = bls12_377_Fq2(bls12_377_Fq("92949345220277864758624960506473182677953048909283248980960104381795901929519566951595905490535835115111760994353"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[2] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410946"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[3] = bls12_377_Fq2(bls12_377_Fq("216465761340224619389371505802605247630151569547285782856803747159100223055385581585702401816380679166954762214499"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[4] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[5] = bls12_377_Fq2(bls12_377_Fq("123516416119946754630746545296132064952198520638002533875843642777304321125866014634106496325844844051843001220146"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[6] = bls12_377_Fq2(bls12_377_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458176"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[7] = bls12_377_Fq2(bls12_377_Fq("165715080792691229252027773188420350858440463845631411558924158284924566418821255823372982649037525009328560463824"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[8] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[9] = bls12_377_Fq2(bls12_377_Fq("42198664672744474621281227892288285906241943207628877683080515507620245292955241189266486323192680957485559243678"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[10] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047232"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[11] = bls12_377_Fq2(bls12_377_Fq("135148009893022339379906188398761468584194992116912126664040619889416147222474808140862391813728516072597320238031"),bls12_377_Fq("0"));

    // Choice of short Weierstrass curve and its twist
    // E(Fq): y^2 = x^3 + 1
    bls12_377_coeff_b = bls12_377_Fq("1");
    // We use a type-D twist here, E'(Fq2): y^2 = x^3 + 1/u
    bls12_377_twist = bls12_377_Fq2(bls12_377_Fq("0"), bls12_377_Fq("1"));
    bls12_377_twist_coeff_b = bls12_377_coeff_b * bls12_377_twist.inverse();

    bls12_377_twist_mul_by_b_c0 = bls12_377_coeff_b * bls12_377_Fq2::non_residue;
    bls12_377_twist_mul_by_b_c1 = bls12_377_coeff_b * bls12_377_Fq2::non_residue;
    bls12_377_twist_mul_by_q_X = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410946"),
                                               bls12_377_Fq("0"));
    bls12_377_twist_mul_by_q_Y = bls12_377_Fq2(bls12_377_Fq("216465761340224619389371505802605247630151569547285782856803747159100223055385581585702401816380679166954762214499"),
                                               bls12_377_Fq("0"));

    // Choice of group G1
    // Identities
    bls12_377_G1::G1_zero = bls12_377_G1(bls12_377_Fq::zero(), bls12_377_Fq::one(), bls12_377_Fq::zero());
    bls12_377_G1::G1_one = bls12_377_G1(bls12_377_Fq("81937999373150964239938255573465948239988671502647976594219695644855304257327692006745978603320413799295628339695"), bls12_377_Fq("241266749859715473739788878240585681733927191168601896383759122102112907357779751001206799952863815012735208165030"), bls12_377_Fq::one());

    // Curve coeffs
    bls12_377_G1::coeff_a = bls12_377_Fq::zero();
    bls12_377_G1::coeff_b = bls12_377_coeff_b;

    // Trace of Frobenius
    bls12_377_trace_of_frobenius = bigint_r("9586122913090633730");

    // Cofactor
    bls12_377_G1::h = bigint<bls12_377_G1::h_limbs>("30631250834960419227450344600217059328");

    // G1 fast subgroup check:  0 == [c0]P + [c1]sigma(P)
    bls12_377_g1_endomorphism_beta = bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945");
    bls12_377_g1_safe_subgroup_check_c1 = bigint_r("91893752504881257701523279626832445441");

    // G1 proof of subgroup: values used to generate x' s.t. [r]x' = x.
    bls12_377_g1_proof_of_safe_subgroup_w = bigint_r("5285428838741532253824584287042945485047145357130994810877");
    bls12_377_g1_proof_of_safe_subgroup_non_member_x = bls12_377_Fq(
        "55791352246783872404788467909907092509364010229903880203689696498787615734938123558571181995209025075818229621722");
    bls12_377_g1_proof_of_safe_subgroup_non_member_y = bls12_377_Fq(
        "174363855833520138229666723484835348689236585013460554444609730120603741818916846216286948728983932214174344518655");

    // WNAF
    //
    // Note to self (AntoineR): Careful with wNAF as it can lead to SCAs: https://eprint.iacr.org/2019/861.pdf
    // Note to self (AntoineR): The GLV patent (https://patents.google.com/patent/US7110538B2/en) expires in 09/2020.
    // As such, efficient techniques for scalar mult. described in https://cryptosith.org/papers/exppair-20130904.pdf
    // can become of interest then.
    //
    // Below we use the same `wnaf_window_table` as used for other curves
    // TODO: Adjust the `wnaf_window_table` and `fixed_base_exp_window_table`
    bls12_377_G1::wnaf_window_table.resize(0);
    bls12_377_G1::wnaf_window_table.push_back(11);
    bls12_377_G1::wnaf_window_table.push_back(24);
    bls12_377_G1::wnaf_window_table.push_back(60);
    bls12_377_G1::wnaf_window_table.push_back(127);

    // Below we use the same `fixed_base_exp_window_table` as used for other curves
    bls12_377_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    bls12_377_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    bls12_377_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    bls12_377_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    bls12_377_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    bls12_377_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    bls12_377_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    bls12_377_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    bls12_377_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    bls12_377_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    bls12_377_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    bls12_377_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    bls12_377_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    bls12_377_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    bls12_377_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    bls12_377_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    bls12_377_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    bls12_377_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    bls12_377_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    bls12_377_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    bls12_377_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    bls12_377_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    bls12_377_G1::fixed_base_exp_window_table.push_back(34552892);

    // Choice of group G2
    // Identities
    bls12_377_G2::G2_zero = bls12_377_G2(bls12_377_Fq2::zero(), bls12_377_Fq2::one(), bls12_377_Fq2::zero());
    bls12_377_G2::G2_one = bls12_377_G2(bls12_377_Fq2(bls12_377_Fq("111583945774695116443911226257823823434468740249883042837745151039122196680777376765707574547389190084887628324746"), bls12_377_Fq("129066980656703085518157301154335215886082112524378686555873161080604845924984124025594590925548060469686767592854")), bls12_377_Fq2(bls12_377_Fq("168863299724668977183029941347596462608978380503965103341003918678547611204475537878680436662916294540335494194722"), bls12_377_Fq("233892497287475762251335351893618429603672921469864392767514552093535653615809913098097380147379993375817193725968")), bls12_377_Fq2::one());

    // Curve twist coeffs
    bls12_377_G2::coeff_a = bls12_377_Fq2::zero();
    bls12_377_G2::coeff_b = bls12_377_twist_coeff_b;

    // Cofactor
    bls12_377_G2::h = bigint<bls12_377_G2::h_limbs>("7923214915284317143930293550643874566881017850177945424769256759165301436616933228209277966774092486467289478618404761412630691835764674559376407658497");

    // Untwist-Frobenius-Twist coefficients
    bls12_377_Fq12 untwist_frobenius_twist_w = bls12_377_Fq12(bls12_377_Fq6::zero(), bls12_377_Fq6::one());
    bls12_377_g2_untwist_frobenius_twist_v = untwist_frobenius_twist_w * untwist_frobenius_twist_w;
    bls12_377_g2_untwist_frobenius_twist_w_3 = untwist_frobenius_twist_w * bls12_377_g2_untwist_frobenius_twist_v;
    bls12_377_g2_untwist_frobenius_twist_v_inverse = bls12_377_g2_untwist_frobenius_twist_v.inverse();
    bls12_377_g2_untwist_frobenius_twist_w_3_inverse = bls12_377_g2_untwist_frobenius_twist_w_3.inverse();

    // Fast cofactor multiplication coefficients
    bls12_377_g2_mul_by_cofactor_h2_0 = bigint_r("293634935485640680722085584138834120318524213360527933441");
    bls12_377_g2_mul_by_cofactor_h2_1 = bigint_r("30631250834960419227450344600217059328");

    // G2 wNAF window table
    bls12_377_G2::wnaf_window_table.resize(0);
    bls12_377_G2::wnaf_window_table.push_back(5);
    bls12_377_G2::wnaf_window_table.push_back(15);
    bls12_377_G2::wnaf_window_table.push_back(39);
    bls12_377_G2::wnaf_window_table.push_back(109);

    // G2 fixed-base exponentiation table
    bls12_377_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    bls12_377_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    bls12_377_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    bls12_377_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    bls12_377_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    bls12_377_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    bls12_377_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    bls12_377_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    bls12_377_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    bls12_377_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    bls12_377_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    bls12_377_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    bls12_377_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    bls12_377_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    bls12_377_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    bls12_377_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    bls12_377_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    bls12_377_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    bls12_377_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    bls12_377_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    bls12_377_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    bls12_377_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    bls12_377_G2::fixed_base_exp_window_table.push_back(31673815);

    // Pairing parameters
    // sage: u = 9586122913090633729
    // sage: ceil(log(u, 2)) # = 64
    // sage: bin(u) # =
    //   '0b1000010100001000110000000000000000000000000000000000000000000001'
    // The Hamming weight of u is: HW(u) = 7, where
    //   u = 2**63 + 2**58 + 2**56 + 2**51 + 2**47 + 2**46 + 1
    // Based on the power-2 decomposition of u, we should have 63 doubling
    // steps and 7 addition steps in the Miller Loop.
    bls12_377_ate_loop_count = bigint_q("9586122913090633729");
    bls12_377_ate_is_loop_count_neg = false;
    // k (embedding degree) = 12
    // bls12_377_final_exponent = (q^12 - 1) / r
    bls12_377_final_exponent = bigint<12*bls12_377_q_limbs>("10623521018019860488254031663707568428798032905123811199571213965079129114663661236359849629341526275899063345613340067081670062620727617884137487754739150147491204559514205186492385590272208934467461444944652711005169371168250068790820776124772095630237102189827733019989835063334551453893534663070786533932633573962932272563471643288531959637300817070265537429506484880990981069041269405383502889677357082012807298529931118124428569059822346289745077401570134157444973271520981774047146918354408632568723153146248333028827919406785654402107153546667815607201488590832478225403444136409349877481268154817904541340614173261949772403060924324366861723245182619859389254985008236007465814273361497134138868945580557938161335670207544906643574043606819537336472235809927599628123275314288006170804044560238676463931639339711913111080974582593228138704154320599775683095604041309000197025419968125718018311805959315220036948621879242495199408833915486421612374480018459896018440926235261824654956932384859260479372776022979736734221629097297890154692194441528462770218811795624471108972377573690833913231260547835550851256817740247389770320334698430697237343583761719223414894063451411431859122738488311580005412765070251810159991897110936324943232526870280724876946523218213525646968094720");
    bls12_377_final_exponent_z = bigint_q("9586122913090633729");
    bls12_377_final_exponent_is_z_neg = false;
}

} // libff
