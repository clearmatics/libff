#include <libff/algebra/curves/bls12_381/bls12_381_g1.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_g2.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_init.hpp>

namespace libff
{

bigint<bls12_381_r_limbs> bls12_381_modulus_r;
bigint<bls12_381_q_limbs> bls12_381_modulus_q;

bls12_381_Fq bls12_381_coeff_b;
bigint<bls12_381_r_limbs> bls12_381_trace_of_frobenius;
bls12_381_Fq2 bls12_381_twist;
bls12_381_Fq2 bls12_381_twist_coeff_b;
bls12_381_Fq bls12_381_twist_mul_by_b_c0;
bls12_381_Fq bls12_381_twist_mul_by_b_c1;
bls12_381_Fq2 bls12_381_twist_mul_by_q_X;
bls12_381_Fq2 bls12_381_twist_mul_by_q_Y;

bigint<bls12_381_r_limbs> bls12_381_g1_safe_subgroup_check_c1;

// Coefficients for G2 untwist-frobenius-twist
bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_v;
bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_w_3;
bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_v_inverse;
bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_w_3_inverse;

bigint<bls12_381_q_limbs> bls12_381_ate_loop_count;
bool bls12_381_ate_is_loop_count_neg;
bigint<12 * bls12_381_q_limbs> bls12_381_final_exponent;
bigint<bls12_381_q_limbs> bls12_381_final_exponent_z;
bool bls12_381_final_exponent_is_z_neg;

void init_bls12_381_params()
{
    using bigint_r = bigint<bls12_381_r_limbs>;
    using bigint_q = bigint<bls12_381_q_limbs>;

    assert(
        sizeof(mp_limb_t) == 8 ||
        sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    bls12_381_modulus_r = bigint_r("5243587517512619047944774050818596583769055"
                                   "2500527637822603658699938581184513");
    assert(bls12_381_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8) {
        bls12_381_Fr::Rsquared =
            bigint_r("329490647479426544212979752063071073927857568219980068178"
                     "8903916070560242797");
        bls12_381_Fr::Rcubed =
            bigint_r("498292539885403193545507422492760844601274463553159150895"
                     "27227471280320770991");
        bls12_381_Fr::inv = 0xfffffffeffffffff; // (-1/modulus) mod W
    }
    if (sizeof(mp_limb_t) == 4) {
        bls12_381_Fr::Rsquared =
            bigint_r("329490647479426544212979752063071073927857568219980068178"
                     "8903916070560242797");
        bls12_381_Fr::Rcubed =
            bigint_r("498292539885403193545507422492760844601274463553159150895"
                     "27227471280320770991");
        bls12_381_Fr::inv = 0xffffffff;
    }
    bls12_381_Fr::num_bits = 255;
    bls12_381_Fr::euler = bigint_r("2621793758756309523972387025409298291884527"
                                   "6250263818911301829349969290592256");
    bls12_381_Fr::s = 32; // 2-adic order of modulus-1
    bls12_381_Fr::t = bigint_r("12208678567578594777604504606729831043093128246"
                               "378069236549469339647"); //(modulus-1)/2^s
    bls12_381_Fr::t_minus_1_over_2 = bigint_r(
        "6104339283789297388802252303364915521546564123189034618274734669823");
    bls12_381_Fr::multiplicative_generator = bls12_381_Fr("7");
    bls12_381_Fr::root_of_unity =
        bls12_381_Fr("102382273577394958236510305758492320625588601802844775411"
                     "89508159991286009131");
    bls12_381_Fr::nqr = bls12_381_Fr("5");
    bls12_381_Fr::nqr_to_t =
        bls12_381_Fr("937917089079007706106976984802249742464848817460758522850"
                     "752807661925904159");
    bls12_381_Fr::static_init();

    /* parameters for base field Fq */
    bls12_381_modulus_q =
        bigint_q("4002409555221667393417789825735904156556882819939007885332058"
                 "136124031650490837864442687629129015664037894272559787");
    assert(bls12_381_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8) {
        bls12_381_Fq::Rsquared = bigint_q(
            "270826391065473017479378762632817651183645519716631767700615429398"
            "2164122222515399004018013397331347120527951271750"); // k=6
        bls12_381_Fq::Rcubed = bigint_q(
            "163906754277462589423671657554808490593875383721159409588363701458"
            "2201460755008380976950835174037649440777609978336");

        bls12_381_Fq::inv = 0x89f3fffcfffcfffd;
    }
    if (sizeof(mp_limb_t) == 4) {
        bls12_381_Fq::Rsquared = bigint_q(
            "270826391065473017479378762632817651183645519716631767700615429398"
            "2164122222515399004018013397331347120527951271750");
        bls12_381_Fq::Rcubed = bigint_q(
            "163906754277462589423671657554808490593875383721159409588363701458"
            "2201460755008380976950835174037649440777609978336");
        bls12_381_Fq::inv = 0xfffcfffd;
    }
    bls12_381_Fq::num_bits = 381;
    bls12_381_Fq::euler =
        bigint_q("2001204777610833696708894912867952078278441409969503942666029"
                 "068062015825245418932221343814564507832018947136279893");
    bls12_381_Fq::s = 1;
    bls12_381_Fq::t =
        bigint_q("2001204777610833696708894912867952078278441409969503942666029"
                 "068062015825245418932221343814564507832018947136279893");
    bls12_381_Fq::t_minus_1_over_2 =
        bigint_q("1000602388805416848354447456433976039139220704984751971333014"
                 "534031007912622709466110671907282253916009473568139946");
    bls12_381_Fq::multiplicative_generator = bls12_381_Fq("2");
    bls12_381_Fq::root_of_unity = bls12_381_Fq(
        "4002409555221667393417789825735904156556882819939007885332058136124031"
        "650490837864442687629129015664037894272559786");
    bls12_381_Fq::nqr = bls12_381_Fq("2");
    bls12_381_Fq::nqr_to_t = bls12_381_Fq(
        "4002409555221667393417789825735904156556882819939007885332058136124031"
        "650490837864442687629129015664037894272559786");
    bls12_381_Fq::static_init();

    /* parameters for twist field Fq2 */
    bls12_381_Fq2::euler = bigint<2 * bls12_381_q_limbs>(
        "8009641123864852705971874322159486308847560049665276329931192268492988"
        "3742456785717003280396510967149874771927700853652655519422698534529681"
        "0010121051821790554650651713590637900820398474816583070927051183888744"
        "9985712996744742684");
    bls12_381_Fq2::s = 3;
    bls12_381_Fq2::t = bigint<2 * bls12_381_q_limbs>(
        "2002410280966213176492968580539871577211890012416319082482798067123247"
        "0935614196429250820099127741787468692981925213413163879855674633632420"
        "2502530262955447638662662928397659475205099618704145767731762795972186"
        "2496428249186185671");
    bls12_381_Fq2::t_minus_1_over_2 = bigint<2 * bls12_381_q_limbs>(
        "1001205140483106588246484290269935788605945006208159541241399033561623"
        "5467807098214625410049563870893734346490962606706581939927837316816210"
        "1251265131477723819331331464198829737602549809352072883865881397986093"
        "1248214124593092835");
    bls12_381_Fq2::non_residue = bls12_381_Fq(
        "4002409555221667393417789825735904156556882819939007885332058136124031"
        "650490837864442687629129015664037894272559786");
    bls12_381_Fq2::nqr =
        bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("1")); // u+1
    bls12_381_Fq2::nqr_to_t = bls12_381_Fq2(
        bls12_381_Fq(
            "102873214623510634997532447921579527738483993692975789615564311803"
            "2610843298655225875571310552543014690878354869257"),
        bls12_381_Fq(
            "297367740898656104344246534652010887917204288300924998917641501809"
            "1420807192182638567116318576472649347015917690530"));
    bls12_381_Fq2::Frobenius_coeffs_c1[0] = bls12_381_Fq("1");
    bls12_381_Fq2::Frobenius_coeffs_c1[1] = bls12_381_Fq(
        "4002409555221667393417789825735904156556882819939007885332058136124031"
        "650490837864442687629129015664037894272559786");
    bls12_381_Fq2::static_init();

    /* parameters for Fq6 */

    bls12_381_Fq6::euler = bigint<6 * bls12_381_q_limbs>(
        "2055413310034685917547178203792332860309200402936847589812920000288065"
        "0925292078031236336437155417867098264344518387365893430771451315477306"
        "7723923870911995741881272581302001907244790574962263945847542833396860"
        "73743682966778056496"
        "2327945395707855461784923849488018385374868097169740055671857527378364"
        "8469851242261268044817816320342666827076722752447529192782544353128746"
        "1471193084577848300836833318729082346882823602164341569076593462295099"
        "0371613607731757827"
        "4075669149520898024347473697702653612215721050521068924301068068177428"
        "7151859717713146107915044671570816889418683602912643322766216203471482"
        "2884004062053629214182533388992931530312083763262100940571236423950189"
        "3128509197213249204");
    bls12_381_Fq6::s = 3;
    bls12_381_Fq6::t = bigint<6 * bls12_381_q_limbs>(
        "5138533275086714793867945509480832150773001007342118974532300000720162"
        "7313230195078090841092888544667745660861295968414733576928628288693266"
        "9309809677279989354703181453255004768111976437405659864618857083492151"
        "8435920741694514124"
        "0581986348926963865446230962372004596343717024292435013917964381844591"
        "2117462810565317011204454080085666706769180688111882298195636088282186"
        "5367798271144462075209208329682270586720705900541085392269148365573774"
        "7592903401932939456"
        "8518917287380224506086868424425663403053930262630267231075267017044357"
        "1787964929428286526978761167892704222354670900728160830691554050867870"
        "5721001015513407303545633347248232882578020940815525235142809105987547"
        "3282127299303312301");
    bls12_381_Fq6::t_minus_1_over_2 = bigint<6 * bls12_381_q_limbs>(
        "2569266637543357396933972754740416075386500503671059487266150000360081"
        "3656615097539045420546444272333872830430647984207366788464314144346633"
        "4654904838639994677351590726627502384055988218702829932309428541746075"
        "9217960370847257062"
        "0290993174463481932723115481186002298171858512146217506958982190922295"
        "6058731405282658505602227040042833353384590344055941149097818044141093"
        "2683899135572231037604604164841135293360352950270542696134574182786887"
        "3796451700966469728"
        "4259458643690112253043434212212831701526965131315133615537633508522178"
        "5893982464714143263489380583946352111177335450364080415345777025433935"
        "2860500507756703651772816673624116441289010470407762617571404552993773"
        "6641063649651656150");
    bls12_381_Fq6::non_residue =
        bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("1"));
    bls12_381_Fq6::nqr = bls12_381_Fq6(
        bls12_381_Fq2::one(), bls12_381_Fq2::one(), bls12_381_Fq2::zero());
    bls12_381_Fq temp_Fq6 = bls12_381_Fq(
        "2973677408986561043442465346520108879172042883009249989176415018091420"
        "807192182638567116318576472649347015917690530");
    bls12_381_Fq6::nqr_to_t = bls12_381_Fq6(
        bls12_381_Fq2(temp_Fq6, temp_Fq6),
        bls12_381_Fq2::zero(),
        bls12_381_Fq2::zero());
    bls12_381_Fq6::Frobenius_coeffs_c1[0] =
        bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c1[1] = bls12_381_Fq2(
        bls12_381_Fq("0"),
        bls12_381_Fq(
            "400240955522166739262431043500668864393550311830558643827117139584"
            "2971157480381377015405980053539358417135540939436"));
    bls12_381_Fq6::Frobenius_coeffs_c1[2] = bls12_381_Fq2(
        bls12_381_Fq("793479390729215512621379701633421447060886740281060493010"
                     "456487427281649075476305620758731620350"),
        bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c1[3] =
        bls12_381_Fq2(bls12_381_Fq("0"), bls12_381_Fq("1"));
    bls12_381_Fq6::Frobenius_coeffs_c1[4] = bls12_381_Fq2(
        bls12_381_Fq(
            "400240955522166739262431043500668864393550311830558643827117139584"
            "2971157480381377015405980053539358417135540939436"),
        bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c1[5] = bls12_381_Fq2(
        bls12_381_Fq("0"),
        bls12_381_Fq("793479390729215512621379701633421447060886740281060493010"
                     "456487427281649075476305620758731620350"));
    bls12_381_Fq6::Frobenius_coeffs_c2[0] =
        bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[1] = bls12_381_Fq2(
        bls12_381_Fq(
            "400240955522166739262431043500668864393550311830558643827117139584"
            "2971157480381377015405980053539358417135540939437"),
        bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[2] = bls12_381_Fq2(
        bls12_381_Fq(
            "400240955522166739262431043500668864393550311830558643827117139584"
            "2971157480381377015405980053539358417135540939436"),
        bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[3] = bls12_381_Fq2(
        bls12_381_Fq(
            "400240955522166739341778982573590415655688281993900788533205813612"
            "4031650490837864442687629129015664037894272559786"),
        bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[4] = bls12_381_Fq2(
        bls12_381_Fq("793479390729215512621379701633421447060886740281060493010"
                     "456487427281649075476305620758731620350"),
        bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[5] = bls12_381_Fq2(
        bls12_381_Fq("793479390729215512621379701633421447060886740281060493010"
                     "456487427281649075476305620758731620351"),
        bls12_381_Fq("0"));

    /* parameters for Fq12 */

    bls12_381_Fq12::euler = bigint<12 * bls12_381_q_limbs>(
        "8449447750135487786386536757818793037762212639342076082205056847631733"
        "7611486762037414270042798116043101246277209976512094361633120467161014"
        "0655544453405655390597288113542772613560509732197200272091941904878386"
        "92659063939320510475"
        "6684686564553744540192748080653129821509884735609423581664936489210967"
        "0440467529068420483514604958227574975684523172972030890413264728910407"
        "5914291222737139398079433126860286068573844453009747181852522430304241"
        "3280382576694575299"
        "4243167774120254336745536274773695755236419654903551716949533916736161"
        "3569172920059467104043034274069321829696761906239463005630083855610676"
        "3024021975665415739763789028633084871968346012196283360672552982454440"
        "6779592032379877481"
        "1619354467715920346231349235866829244117661464414438262378914602649540"
        "0686951521638441608423046735329190594695892803338108752058997137128797"
        "3778883772765769448885014782339635406211030364946185240866112535299846"
        "0333970628572435583"
        "1757875201144654542881444863180645429521461839597984159381168059700378"
        "1998317840382782546961827784098634063653000816109224177006797686172283"
        "6524688010589950951561114449900951002753319232147895757766930045754791"
        "8189673190852298371"
        "5647603567520167155789018111526068121587518840305926270907547770597694"
        "1239754628313695901797410453781873141027101200728931074856266540267317"
        "2430117685618695762741237977291546454709306808974723155022101617076402"
        "8215685592439765640");
    bls12_381_Fq12::s = 4;
    bls12_381_Fq12::t = bigint<12 * bls12_381_q_limbs>(
        "1056180968766935973298317094727349129720276579917759510275632105953966"
        "7201435845254676783755349764505387655784651247064011795204140058395126"
        "7581943056675706923824661014192846576695063716524650034011492738109798"
        "36582382992415063809"
        "4585585820569218067524093510081641227688735591951177947708117061151370"
        "8805058441133552560439325619778446871960565396621503861301658091113800"
        "9489286402842142424759929140857535758571730556626218397731565303788030"
        "1660047822086821912"
        "4280395971765031792093192034346711969404552456862943964618691739592020"
        "1696146615007433388005379284258665228712095238279932875703760481951334"
        "5378002746958176967470473628579135608996043251524535420084069122806805"
        "0847449004047484685"
        "1452419308464490043278918654483353655514707683051804782797364325331192"
        "5085868940204805201052880841916148824336986600417263594007374642141099"
        "6722360471595721181110626847792454425776378795618273155108264066912480"
        "7541746328571554447"
        "8969734400143081817860180607897580678690182729949748019922646007462547"
        "2749789730047847818370228473012329257956625102013653022125849710771535"
        "4565586001323743868945139306237618875344164904018486969720866255719348"
        "9773709148856537296"
        "4455950445940020894473627263940758515198439855038240783863443471324711"
        "7654969328539211987724676306722734142628387650091116384357033317533414"
        "6553764710702336970342654747161443306838663351121840394377762702134550"
        "3526960699054970705");
    bls12_381_Fq12::t_minus_1_over_2 = bigint<12 * bls12_381_q_limbs>(
        "5280904843834679866491585473636745648601382899588797551378160529769833"
        "6007179226273383918776748822526938278923256235320058976020700291975633"
        "7909715283378534619123305070964232883475318582623250170057463690548991"
        "8291191496207531904"
        "7292792910284609033762046755040820613844367795975588973854058530575685"
        "4402529220566776280219662809889223435980282698310751930650829045556900"
        "4744643201421071212379964570428767879285865278313109198865782651894015"
        "0830023911043410956"
        "2140197985882515896046596017173355984702276228431471982309345869796010"
        "0848073307503716694002689642129332614356047619139966437851880240975667"
        "2689001373479088483735236814289567804498021625762267710042034561403402"
        "5423724502023742342"
        "5726209654232245021639459327241676827757353841525902391398682162665596"
        "2542934470102402600526440420958074412168493300208631797003687321070549"
        "8361180235797860590555313423896227212888189397809136577554132033456240"
        "3770873164285777223"
        "9484867200071540908930090303948790339345091364974874009961323003731273"
        "6374894865023923909185114236506164628978312551006826511062924855385767"
        "7282793000661871934472569653118809437672082452009243484860433127859674"
        "4886854574428268648"
        "2227975222970010447236813631970379257599219927519120391931721735662355"
        "8827484664269605993862338153361367071314193825045558192178516658766707"
        "3276882355351168485171327373580721653419331675560920197188881351067275"
        "1763480349527485352");
    bls12_381_Fq12::non_residue =
        bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("1"));
    bls12_381_Fq12::nqr =
        bls12_381_Fq12(bls12_381_Fq6::zero(), bls12_381_Fq6::one());
    bls12_381_Fq temp_Fq12 = bls12_381_Fq(
        "3357996710086603428986649435961018971596863377125478091385687488711898"
        "724126407611022502097010210262797519903698974");
    bls12_381_Fq12::nqr_to_t = bls12_381_Fq12(
        bls12_381_Fq6::zero(),
        bls12_381_Fq6(
            bls12_381_Fq2::zero(),
            bls12_381_Fq2(temp_Fq12, temp_Fq12),
            bls12_381_Fq2::zero()));
    bls12_381_Fq12::Frobenius_coeffs_c1[0] =
        bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[1] = bls12_381_Fq2(
        bls12_381_Fq(
            "385075437003716901195214707605136405715880742097068243867605052261"
            "3628423219637725072182697113062777891589506424760"),
        bls12_381_Fq(
            "151655185184498381465642749684540099398075398968325446656007613510"
            "403227271200139370504932015952886146304766135027"));
    bls12_381_Fq12::Frobenius_coeffs_c1[2] = bls12_381_Fq2(
        bls12_381_Fq("793479390729215512621379701633421447060886740281060493010"
                     "456487427281649075476305620758731620351"),
        bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[3] = bls12_381_Fq2(
        bls12_381_Fq(
            "297367740898656104344246534652010887917204288300924998917641501809"
            "1420807192182638567116318576472649347015917690530"),
        bls12_381_Fq(
            "102873214623510634997532447921579527738483993692975789615564311803"
            "2610843298655225875571310552543014690878354869257"));
    bls12_381_Fq12::Frobenius_coeffs_c1[4] = bls12_381_Fq2(
        bls12_381_Fq("793479390729215512621379701633421447060886740281060493010"
                     "456487427281649075476305620758731620350"),
        bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[5] = bls12_381_Fq2(
        bls12_381_Fq(
            "312533259417105942490810809620464897857011828197757543583242263160"
            "1824034463382777937621250592425535493320683825557"),
        bls12_381_Fq(
            "877076961050607968509681729531255177986764537961432449499635504522"
            "207616027455086505066378536590128544573588734230"));
    bls12_381_Fq12::Frobenius_coeffs_c1[6] = bls12_381_Fq2(
        bls12_381_Fq(
            "400240955522166739341778982573590415655688281993900788533205813612"
            "4031650490837864442687629129015664037894272559786"),
        bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[7] = bls12_381_Fq2(
        bls12_381_Fq(
            "151655185184498381465642749684540099398075398968325446656007613510"
            "403227271200139370504932015952886146304766135027"),
        bls12_381_Fq(
            "385075437003716901195214707605136405715880742097068243867605052261"
            "3628423219637725072182697113062777891589506424760"));
    bls12_381_Fq12::Frobenius_coeffs_c1[8] = bls12_381_Fq2(
        bls12_381_Fq(
            "400240955522166739262431043500668864393550311830558643827117139584"
            "2971157480381377015405980053539358417135540939436"),
        bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[9] = bls12_381_Fq2(
        bls12_381_Fq(
            "102873214623510634997532447921579527738483993692975789615564311803"
            "2610843298655225875571310552543014690878354869257"),
        bls12_381_Fq(
            "297367740898656104344246534652010887917204288300924998917641501809"
            "1420807192182638567116318576472649347015917690530"));
    bls12_381_Fq12::Frobenius_coeffs_c1[10] = bls12_381_Fq2(
        bls12_381_Fq(
            "400240955522166739262431043500668864393550311830558643827117139584"
            "2971157480381377015405980053539358417135540939437"),
        bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[11] = bls12_381_Fq2(
        bls12_381_Fq(
            "877076961050607968509681729531255177986764537961432449499635504522"
            "207616027455086505066378536590128544573588734230"),
        bls12_381_Fq(
            "312533259417105942490810809620464897857011828197757543583242263160"
            "1824034463382777937621250592425535493320683825557"));

    // Choice of short Weierstrass curve and its twist
    // E(Fq): y^2 = x^3 + 4

    bls12_381_coeff_b = bls12_381_Fq("4");
    bls12_381_twist = bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("1"));
    bls12_381_twist_coeff_b = bls12_381_coeff_b * bls12_381_twist;
    bls12_381_twist_mul_by_b_c0 =
        bls12_381_coeff_b * bls12_381_Fq2::non_residue;
    bls12_381_twist_mul_by_b_c1 =
        bls12_381_coeff_b * bls12_381_Fq2::non_residue;
    bls12_381_twist_mul_by_q_X = bls12_381_Fq2(
        bls12_381_Fq("0"),
        bls12_381_Fq(
            "400240955522166739262431043500668864393550311830558643827117139584"
            "2971157480381377015405980053539358417135540939437"));
    bls12_381_twist_mul_by_q_Y = bls12_381_Fq2(
        bls12_381_Fq(
            "297367740898656104344246534652010887917204288300924998917641501809"
            "1420807192182638567116318576472649347015917690530"),
        bls12_381_Fq(
            "102873214623510634997532447921579527738483993692975789615564311803"
            "2610843298655225875571310552543014690878354869257"));

    /* choice of group G1 */
    bls12_381_G1::G1_zero = bls12_381_G1(
        bls12_381_Fq::zero(), bls12_381_Fq::one(), bls12_381_Fq::zero());
    bls12_381_G1::G1_one = bls12_381_G1(
        bls12_381_Fq(
            "368541675371338701678108831518307775796162079578254640989457837868"
            "8607592378376318836054947676345821548104185464507"),
        bls12_381_Fq(
            "133950654494447647302047137994192122158493387593834962042654373641"
            "6511423956333506472724655353366534992391756441569"),
        bls12_381_Fq::one());

    // Curve coeffs
    bls12_381_G1::coeff_a = bls12_381_Fq::zero();
    bls12_381_G1::coeff_b = bls12_381_coeff_b;

    // Cofactor
    bls12_381_G1::h =
        bigint<bls12_381_G1::h_limbs>("76329603384216526031706109802092473003");

    // TODO: wNAF window table
    bls12_381_G1::wnaf_window_table.resize(0);
    bls12_381_G1::wnaf_window_table.push_back(11);
    bls12_381_G1::wnaf_window_table.push_back(24);
    bls12_381_G1::wnaf_window_table.push_back(60);
    bls12_381_G1::wnaf_window_table.push_back(127);

    // TODO: fixed-base exponentiation table
    bls12_381_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    bls12_381_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    bls12_381_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    bls12_381_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    bls12_381_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    bls12_381_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    bls12_381_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    bls12_381_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    bls12_381_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    bls12_381_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    bls12_381_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    bls12_381_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    bls12_381_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    bls12_381_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    bls12_381_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    bls12_381_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    bls12_381_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    bls12_381_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    bls12_381_G1::fixed_base_exp_window_table.push_back(34552892);

    /* choice of group G2 */
    bls12_381_G2::G2_zero = bls12_381_G2(
        bls12_381_Fq2::zero(), bls12_381_Fq2::one(), bls12_381_Fq2::zero());

    // simple G2 generator
    bls12_381_G2::G2_one = bls12_381_G2(
        bls12_381_Fq2(
            bls12_381_Fq(
                "35270106958746661818713911601106014489002995279277524021990864"
                "4239793785735715026873347600343865175952761926303160"),
            bls12_381_Fq(
                "30591443442442137099712598147537816369864703254766475586593732"
                "06291635324768958432433509563104347017837885763365758")),
        bls12_381_Fq2(
            bls12_381_Fq(
                "19851506022872919355680545211771716383008689782156557308593786"
                "65066344726373823718423869104263333984641494340347905"),
            bls12_381_Fq(
                "92755366549233245574720196577603788075774019345359297002502797"
                "8793976877002675564980949289727957565575433344219582")),
        bls12_381_Fq2::one());

    // Curve twist coeffs
    bls12_381_G2::coeff_a = bls12_381_Fq2::zero();
    bls12_381_G2::coeff_b = bls12_381_twist_coeff_b;

    // Cofactor
    bls12_381_G2::h = bigint<bls12_381_G2::h_limbs>(
        "3055023339312683442009997531931215042144660192541881426676640329822676"
        "0418297188402650742735925997784783227283904161666128580382337837209635"
        "5777062779109");

    // TODO: wNAF window table
    bls12_381_G2::wnaf_window_table.resize(0);
    bls12_381_G2::wnaf_window_table.push_back(5);
    bls12_381_G2::wnaf_window_table.push_back(15);
    bls12_381_G2::wnaf_window_table.push_back(39);
    bls12_381_G2::wnaf_window_table.push_back(109);

    // TODO: fixed-base exponentiation table
    bls12_381_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    bls12_381_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    bls12_381_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    bls12_381_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    bls12_381_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    bls12_381_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    bls12_381_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    bls12_381_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    bls12_381_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    bls12_381_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    bls12_381_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    bls12_381_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    bls12_381_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    bls12_381_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    bls12_381_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    bls12_381_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    bls12_381_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    bls12_381_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    bls12_381_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    bls12_381_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    bls12_381_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    bls12_381_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    bls12_381_G2::fixed_base_exp_window_table.push_back(31673815);

    /* pairing parameters */

    bls12_381_ate_loop_count =
        bigint<bls12_381_q_limbs>("15132376222941642752");
    bls12_381_ate_is_loop_count_neg = true;
    bls12_381_final_exponent = bigint<12 * bls12_381_q_limbs>(
        "3222773615169341404628915645865101399083799695148284942183666880252886"
        "6104110468279499868049758000889997324981410444769277898820837677957381"
        "9485263026159588510513834876303014016798809919343532899164848730280942"
        "6099566709175656181158672873996232868132703579017315101881499343633603"
        "8161450133408682544227192007936328995451056537537844370437299488140679"
        "7882676971082200626541916413184642520269678897559532260949334760604962"
        "0863488981189822488426343796375986654688177690758785554937522144927901"
        "2278585020295757520017608420442275148595733646547232481098283363849090"
        "4279282696134323072515220044451592646885410572234451732790590013479358"
        "3438412200741748482217220170835978720176385141031741227848439255783704"
        "3084352295960009567628572373704943834654475316891297497679152853527631"
        "7256904336520179281145394686565050419250614107803233314658825463117900"
        "2507011991815292059423631593257659918194339143039088604607205814082013"
        "7316404777379482541101192230582006561112154456180841405530221205747139"
        "5719432072209245600258134364584636810093520285711072578721435517884103"
        "5264838327332898024261573015427444767400084947803633543051169788056206"
        "7146707140071135883955337534072489973546048014459978201490658654381329"
        "2157922220645089192130209334926661588737007768565838519456601560804957"
        "985667880395221049249803753582637708560");
    bls12_381_final_exponent_z =
        bigint<bls12_381_q_limbs>("15132376222941642752");
    bls12_381_final_exponent_is_z_neg = true;
}

} // namespace libff
