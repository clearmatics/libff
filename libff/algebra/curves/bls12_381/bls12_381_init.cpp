#include <libff/algebra/curves/bls12_381/bls12_381_g1.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_g2.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_init.hpp>

namespace libff
{

bigint<bls12_381_r_limbs> bls12_381_modulus_r; // (VV)

bls12_381_Fq bls12_381_coeff_b;
bigint<bls12_381_r_limbs> bls12_381_trace_of_frobenius; // from bls12_377 (VV)
bls12_381_Fq2 bls12_381_twist;
bls12_381_Fq2 bls12_381_twist_coeff_b;
bls12_381_Fq bls12_381_twist_mul_by_b_c0;
bls12_381_Fq bls12_381_twist_mul_by_b_c1;
bls12_381_Fq2 bls12_381_twist_mul_by_q_X;
bls12_381_Fq2 bls12_381_twist_mul_by_q_Y;

// from bls12_377 (VV) 
bls12_381_Fq bls12_381_g1_endomorphism_beta;
bigint<bls12_381_r_limbs> bls12_381_g1_safe_subgroup_check_c1;


// Coefficients for G2 untwist-frobenius-twist (from bls12_377 VV)
bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_v;
bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_w_3;
bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_v_inverse;
bls12_381_Fq12 bls12_381_g2_untwist_frobenius_twist_w_3_inverse;

  
bigint<bls12_381_q_limbs> bls12_381_ate_loop_count;
bool bls12_381_ate_is_loop_count_neg;
bigint<12*bls12_381_q_limbs> bls12_381_final_exponent;
bigint<bls12_381_q_limbs> bls12_381_final_exponent_z;
bool bls12_381_final_exponent_is_z_neg;

void init_bls12_381_params()
{
    //    init_bls12_381_fields(); (VV)
    using bigint_r = bigint<bls12_381_r_limbs>;
    using bigint_q = bigint<bls12_381_q_limbs>;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    bls12_381_modulus_r = bigint_r("52435875175126190479447740508185965837690552500527637822603658699938581184513");
    assert(bls12_381_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        bls12_381_Fr::Rsquared = bigint_r("3294906474794265442129797520630710739278575682199800681788903916070560242797");
        bls12_381_Fr::Rcubed = bigint_r("49829253988540319354550742249276084460127446355315915089527227471280320770991");
        bls12_381_Fr::inv = 0xfffffffeffffffff; // (-1/modulus) mod W
    }
    if (sizeof(mp_limb_t) == 4)
    {
        bls12_381_Fr::Rsquared = bigint_r("3294906474794265442129797520630710739278575682199800681788903916070560242797");
        bls12_381_Fr::Rcubed = bigint_r("49829253988540319354550742249276084460127446355315915089527227471280320770991");
        bls12_381_Fr::inv = 0xffffffff;
    }
    bls12_381_Fr::num_bits = 255;
    bls12_381_Fr::euler = bigint_r("26217937587563095239723870254092982918845276250263818911301829349969290592256");
    bls12_381_Fr::s = 32; // 2-adic order of modulus-1
    bls12_381_Fr::t = bigint_r("12208678567578594777604504606729831043093128246378069236549469339647"); //(modulus-1)/2^s
    bls12_381_Fr::t_minus_1_over_2 = bigint_r("6104339283789297388802252303364915521546564123189034618274734669823");
    bls12_381_Fr::multiplicative_generator = bls12_381_Fr("7");
    bls12_381_Fr::root_of_unity = bls12_381_Fr("10238227357739495823651030575849232062558860180284477541189508159991286009131");
    bls12_381_Fr::nqr = bls12_381_Fr("5");
    bls12_381_Fr::nqr_to_t = bls12_381_Fr("937917089079007706106976984802249742464848817460758522850752807661925904159");

    /* parameters for base field Fq */
    bls12_381_modulus_q = bigint_q("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787");
    assert(bls12_381_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        bls12_381_Fq::Rsquared = bigint_q("2708263910654730174793787626328176511836455197166317677006154293982164122222515399004018013397331347120527951271750"); // k=6
        bls12_381_Fq::Rcubed = bigint_q("1639067542774625894236716575548084905938753837211594095883637014582201460755008380976950835174037649440777609978336");

        bls12_381_Fq::inv = 0x89f3fffcfffcfffd;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        bls12_381_Fq::Rsquared = bigint_q("2708263910654730174793787626328176511836455197166317677006154293982164122222515399004018013397331347120527951271750");
        bls12_381_Fq::Rcubed = bigint_q("1639067542774625894236716575548084905938753837211594095883637014582201460755008380976950835174037649440777609978336");
        bls12_381_Fq::inv = 0xfffcfffd;
    }
    bls12_381_Fq::num_bits = 381;
    bls12_381_Fq::euler = bigint_q("2001204777610833696708894912867952078278441409969503942666029068062015825245418932221343814564507832018947136279893");
    bls12_381_Fq::s = 1;
    bls12_381_Fq::t = bigint_q("2001204777610833696708894912867952078278441409969503942666029068062015825245418932221343814564507832018947136279893");
    bls12_381_Fq::t_minus_1_over_2 = bigint_q("1000602388805416848354447456433976039139220704984751971333014534031007912622709466110671907282253916009473568139946");
    bls12_381_Fq::multiplicative_generator = bls12_381_Fq("2");
    bls12_381_Fq::root_of_unity = bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786");
    bls12_381_Fq::nqr = bls12_381_Fq("2");
    bls12_381_Fq::nqr_to_t = bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786");

    /* parameters for twist field Fq2 */
    bls12_381_Fq2::euler = bigint<2*bls12_381_q_limbs>("8009641123864852705971874322159486308847560049665276329931192268492988374245678571700328039651096714987477192770085365265551942269853452968100101210518217905546506517135906379008203984748165830709270511838887449985712996744742684");
    bls12_381_Fq2::s = 3;
    bls12_381_Fq2::t = bigint<2*bls12_381_q_limbs>("2002410280966213176492968580539871577211890012416319082482798067123247093561419642925082009912774178746869298192521341316387985567463363242025025302629554476386626629283976594752050996187041457677317627959721862496428249186185671");
    bls12_381_Fq2::t_minus_1_over_2 = bigint<2*bls12_381_q_limbs>("1001205140483106588246484290269935788605945006208159541241399033561623546780709821462541004956387089373434649096260670658193992783731681621012512651314777238193313314641988297376025498093520728838658813979860931248214124593092835");
    bls12_381_Fq2::non_residue = bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786");
    bls12_381_Fq2::nqr = bls12_381_Fq2(bls12_381_Fq("1"),bls12_381_Fq("1")); // u+1
    bls12_381_Fq2::nqr_to_t = bls12_381_Fq2(bls12_381_Fq("1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"),bls12_381_Fq("2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530"));
    bls12_381_Fq2::Frobenius_coeffs_c1[0] = bls12_381_Fq("1");
    bls12_381_Fq2::Frobenius_coeffs_c1[1] = bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786");

    /* parameters for Fq6 */

    bls12_381_Fq6::euler = bigint<6*bls12_381_q_limbs>("20554133100346859175471782037923328603092004029368475898129200002880650925292078031236336437155417867098264344518387365893430771451315477306772392387091199574188127258130200190724479057496226394584754283339686073743682966778056496"
                                                       "2327945395707855461784923849488018385374868097169740055671857527378364846985124226126804481781632034266682707672275244752919278254435312874614711930845778483008368333187290823468828236021643415690765934622950990371613607731757827"
                                                       "4075669149520898024347473697702653612215721050521068924301068068177428715185971771314610791504467157081688941868360291264332276621620347148228840040620536292141825333889929315303120837632621009405712364239501893128509197213249204");
    bls12_381_Fq6::s = 3;
    bls12_381_Fq6::t = bigint<6*bls12_381_q_limbs>("5138533275086714793867945509480832150773001007342118974532300000720162731323019507809084109288854466774566086129596841473357692862828869326693098096772799893547031814532550047681119764374056598646188570834921518435920741694514124"
                                                   "0581986348926963865446230962372004596343717024292435013917964381844591211746281056531701120445408008566670676918068811188229819563608828218653677982711444620752092083296822705867207059005410853922691483655737747592903401932939456"
                                                   "8518917287380224506086868424425663403053930262630267231075267017044357178796492942828652697876116789270422235467090072816083069155405086787057210010155134073035456333472482328825780209408155252351428091059875473282127299303312301");
    bls12_381_Fq6::t_minus_1_over_2 = bigint<6*bls12_381_q_limbs>("2569266637543357396933972754740416075386500503671059487266150000360081365661509753904542054644427233387283043064798420736678846431414434663346549048386399946773515907266275023840559882187028299323094285417460759217960370847257062"
                                                                  "0290993174463481932723115481186002298171858512146217506958982190922295605873140528265850560222704004283335338459034405594114909781804414109326838991355722310376046041648411352933603529502705426961345741827868873796451700966469728"
                                                                  "4259458643690112253043434212212831701526965131315133615537633508522178589398246471414326348938058394635211117733545036408041534577702543393528605005077567036517728166736241164412890104704077626175714045529937736641063649651656150");
    bls12_381_Fq6::non_residue = bls12_381_Fq2(bls12_381_Fq("1"),bls12_381_Fq("1"));
    bls12_381_Fq6::nqr = bls12_381_Fq6(bls12_381_Fq2::one(),bls12_381_Fq2::one(),bls12_381_Fq2::zero());
    bls12_381_Fq temp_Fq6 = bls12_381_Fq("2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530");
    bls12_381_Fq6::nqr_to_t = bls12_381_Fq6(bls12_381_Fq2(temp_Fq6,temp_Fq6),bls12_381_Fq2::zero(),bls12_381_Fq2::zero());
    bls12_381_Fq6::Frobenius_coeffs_c1[0] = bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c1[1] = bls12_381_Fq2(bls12_381_Fq("0"), bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436"));
    bls12_381_Fq6::Frobenius_coeffs_c1[2] = bls12_381_Fq2(bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c1[3] = bls12_381_Fq2(bls12_381_Fq("0"), bls12_381_Fq("1"));
    bls12_381_Fq6::Frobenius_coeffs_c1[4] = bls12_381_Fq2(bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c1[5] = bls12_381_Fq2(bls12_381_Fq("0"), bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350"));
    bls12_381_Fq6::Frobenius_coeffs_c2[0] = bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[1] = bls12_381_Fq2(bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939437"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[2] = bls12_381_Fq2(bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[3] = bls12_381_Fq2(bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[4] = bls12_381_Fq2(bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[5] = bls12_381_Fq2(bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620351"), bls12_381_Fq("0"));

    /* parameters for Fq12 */

    bls12_381_Fq12::euler = bigint<12*bls12_381_q_limbs>("84494477501354877863865367578187930377622126393420760822050568476317337611486762037414270042798116043101246277209976512094361633120467161014065554445340565539059728811354277261356050973219720027209194190487838692659063939320510475"
                                                         "6684686564553744540192748080653129821509884735609423581664936489210967044046752906842048351460495822757497568452317297203089041326472891040759142912227371393980794331268602860685738444530097471818525224303042413280382576694575299"
                                                         "4243167774120254336745536274773695755236419654903551716949533916736161356917292005946710404303427406932182969676190623946300563008385561067630240219756654157397637890286330848719683460121962833606725529824544406779592032379877481"
                                                         "1619354467715920346231349235866829244117661464414438262378914602649540068695152163844160842304673532919059469589280333810875205899713712879737788837727657694488850147823396354062110303649461852408661125352998460333970628572435583"
                                                         "1757875201144654542881444863180645429521461839597984159381168059700378199831784038278254696182778409863406365300081610922417700679768617228365246880105899509515611144499009510027533192321478957577669300457547918189673190852298371"
                                                         "5647603567520167155789018111526068121587518840305926270907547770597694123975462831369590179741045378187314102710120072893107485626654026731724301176856186957627412379772915464547093068089747231550221016170764028215685592439765640");
    bls12_381_Fq12::s = 4;
    bls12_381_Fq12::t = bigint<12*bls12_381_q_limbs>("10561809687669359732983170947273491297202765799177595102756321059539667201435845254676783755349764505387655784651247064011795204140058395126758194305667570692382466101419284657669506371652465003401149273810979836582382992415063809"
                                                     "4585585820569218067524093510081641227688735591951177947708117061151370880505844113355256043932561977844687196056539662150386130165809111380094892864028421424247599291408575357585717305566262183977315653037880301660047822086821912"
                                                     "4280395971765031792093192034346711969404552456862943964618691739592020169614661500743338800537928425866522871209523827993287570376048195133453780027469581769674704736285791356089960432515245354200840691228068050847449004047484685"
                                                     "1452419308464490043278918654483353655514707683051804782797364325331192508586894020480520105288084191614882433698660041726359400737464214109967223604715957211811106268477924544257763787956182731551082640669124807541746328571554447"
                                                     "8969734400143081817860180607897580678690182729949748019922646007462547274978973004784781837022847301232925795662510201365302212584971077153545655860013237438689451393062376188753441649040184869697208662557193489773709148856537296"
                                                     "4455950445940020894473627263940758515198439855038240783863443471324711765496932853921198772467630672273414262838765009111638435703331753341465537647107023369703426547471614433068386633511218403943777627021345503526960699054970705");
    bls12_381_Fq12::t_minus_1_over_2 = bigint<12*bls12_381_q_limbs>("5280904843834679866491585473636745648601382899588797551378160529769833600717922627338391877674882252693827892325623532005897602070029197563379097152833785346191233050709642328834753185826232501700574636905489918291191496207531904"
                                                                    "7292792910284609033762046755040820613844367795975588973854058530575685440252922056677628021966280988922343598028269831075193065082904555690047446432014210712123799645704287678792858652783131091988657826518940150830023911043410956"
                                                                    "2140197985882515896046596017173355984702276228431471982309345869796010084807330750371669400268964212933261435604761913996643785188024097566726890013734790884837352368142895678044980216257622677100420345614034025423724502023742342"
                                                                    "5726209654232245021639459327241676827757353841525902391398682162665596254293447010240260052644042095807441216849330020863179700368732107054983611802357978605905553134238962272128881893978091365775541320334562403770873164285777223"
                                                                    "9484867200071540908930090303948790339345091364974874009961323003731273637489486502392390918511423650616462897831255100682651106292485538576772827930006618719344725696531188094376720824520092434848604331278596744886854574428268648"
                                                                    "2227975222970010447236813631970379257599219927519120391931721735662355882748466426960599386233815336136707131419382504555819217851665876670732768823553511684851713273735807216534193316755609201971888813510672751763480349527485352");
    bls12_381_Fq12::non_residue = bls12_381_Fq2(bls12_381_Fq("1"),bls12_381_Fq("1"));
    bls12_381_Fq12::nqr = bls12_381_Fq12(bls12_381_Fq6::zero(),bls12_381_Fq6::one());
    bls12_381_Fq temp_Fq12 = bls12_381_Fq("3357996710086603428986649435961018971596863377125478091385687488711898724126407611022502097010210262797519903698974");
    bls12_381_Fq12::nqr_to_t = bls12_381_Fq12(bls12_381_Fq6::zero(),bls12_381_Fq6(bls12_381_Fq2::zero(),bls12_381_Fq2(temp_Fq12,temp_Fq12),bls12_381_Fq2::zero()));
    bls12_381_Fq12::Frobenius_coeffs_c1[0] = bls12_381_Fq2(bls12_381_Fq("1"),bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[1] = bls12_381_Fq2(bls12_381_Fq("3850754370037169011952147076051364057158807420970682438676050522613628423219637725072182697113062777891589506424760"),bls12_381_Fq("151655185184498381465642749684540099398075398968325446656007613510403227271200139370504932015952886146304766135027"));
    bls12_381_Fq12::Frobenius_coeffs_c1[2] = bls12_381_Fq2(bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620351"),bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[3] = bls12_381_Fq2(bls12_381_Fq("2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530"),bls12_381_Fq("1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"));
    bls12_381_Fq12::Frobenius_coeffs_c1[4] = bls12_381_Fq2(bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350"),bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[5] = bls12_381_Fq2(bls12_381_Fq("3125332594171059424908108096204648978570118281977575435832422631601824034463382777937621250592425535493320683825557"),bls12_381_Fq("877076961050607968509681729531255177986764537961432449499635504522207616027455086505066378536590128544573588734230"));
    bls12_381_Fq12::Frobenius_coeffs_c1[6] = bls12_381_Fq2(bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786"),bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[7] = bls12_381_Fq2(bls12_381_Fq("151655185184498381465642749684540099398075398968325446656007613510403227271200139370504932015952886146304766135027"),bls12_381_Fq("3850754370037169011952147076051364057158807420970682438676050522613628423219637725072182697113062777891589506424760"));
    bls12_381_Fq12::Frobenius_coeffs_c1[8] = bls12_381_Fq2(bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436"),bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[9] = bls12_381_Fq2(bls12_381_Fq("1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"),bls12_381_Fq("2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530"));
    bls12_381_Fq12::Frobenius_coeffs_c1[10] = bls12_381_Fq2(bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939437"),bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[11] = bls12_381_Fq2(bls12_381_Fq("877076961050607968509681729531255177986764537961432449499635504522207616027455086505066378536590128544573588734230"),bls12_381_Fq("3125332594171059424908108096204648978570118281977575435832422631601824034463382777937621250592425535493320683825557"));

    /* choice of short Weierstrass curve and its twist */

    bls12_381_coeff_b = bls12_381_Fq("4");
    bls12_381_twist = bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("1"));
    bls12_381_twist_coeff_b = bls12_381_coeff_b * bls12_381_twist;
    bls12_381_twist_mul_by_b_c0 = bls12_381_coeff_b * bls12_381_Fq2::non_residue;
    bls12_381_twist_mul_by_b_c1 = bls12_381_coeff_b * bls12_381_Fq2::non_residue;
    bls12_381_twist_mul_by_q_X = bls12_381_Fq2(bls12_381_Fq("0"),
                                               bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939437"));
    bls12_381_twist_mul_by_q_Y = bls12_381_Fq2(bls12_381_Fq("2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530"),
                                               bls12_381_Fq("1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"));


    /* choice of group G1 */
    bls12_381_G1::G1_zero = bls12_381_G1(bls12_381_Fq::zero(),
                                     bls12_381_Fq::one(),
                                     bls12_381_Fq::zero());
    bls12_381_G1::G1_one = bls12_381_G1(bls12_381_Fq("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507"),
                                    bls12_381_Fq("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569"),
                                    bls12_381_Fq::one());
    // Curve coeffs (from bls12_377 (VV))
    bls12_381_G1::coeff_a = bls12_381_Fq::zero();
    bls12_381_G1::coeff_b = bls12_381_coeff_b;
    
    // Cofactor
    bls12_381_G1::h = bigint<bls12_381_G1::h_limbs>("76329603384216526031706109802092473003");

    // from bls12_377 (VV) !!! change constants for 381
    // G1 fast subgroup check:  0 == [c0]P + [c1]sigma(P)
    bls12_381_g1_endomorphism_beta =
        bls12_381_Fq("809496482649127194085583631406374772648452947207104994781"
                     "37287262712535938301461879813459410945");
    bls12_381_g1_safe_subgroup_check_c1 =
        bigint_r("91893752504881257701523279626832445441");

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
    bls12_381_G2::G2_zero = bls12_381_G2(bls12_381_Fq2::zero(),
                                         bls12_381_Fq2::one(),
                                         bls12_381_Fq2::zero());

    // simple G2 generator
    bls12_381_G2::G2_one = bls12_381_G2(bls12_381_Fq2(bls12_381_Fq("352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160"),
                                                      bls12_381_Fq("3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758")),
                                        bls12_381_Fq2(bls12_381_Fq("1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905"),
                                                      bls12_381_Fq("927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582")),
                                        bls12_381_Fq2::one());
    // Curve twist coeffs (from bls12_377 VV)
    bls12_381_G2::coeff_a = bls12_381_Fq2::zero();
    bls12_381_G2::coeff_b = bls12_381_twist_coeff_b;
    
    // Cofactor
    bls12_381_G2::h = bigint<bls12_381_G2::h_limbs>("305502333931268344200999753193121504214466019254188142667664032982267604182971884026507427359259977847832272839041616661285803823378372096355777062779109");

     // Untwist-Frobenius-Twist coefficients (from bls12_377 VV)
    bls12_381_Fq12 untwist_frobenius_twist_w =
        bls12_381_Fq12(bls12_381_Fq6::zero(), bls12_381_Fq6::one());
    bls12_381_g2_untwist_frobenius_twist_v =
        untwist_frobenius_twist_w * untwist_frobenius_twist_w;
    bls12_381_g2_untwist_frobenius_twist_w_3 =
        untwist_frobenius_twist_w * bls12_381_g2_untwist_frobenius_twist_v;
    bls12_381_g2_untwist_frobenius_twist_v_inverse =
        bls12_381_g2_untwist_frobenius_twist_v.inverse();
    bls12_381_g2_untwist_frobenius_twist_w_3_inverse =
        bls12_381_g2_untwist_frobenius_twist_w_3.inverse();
   

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

    bls12_381_ate_loop_count = bigint<bls12_381_q_limbs>("15132376222941642752");
    bls12_381_ate_is_loop_count_neg = true;
    bls12_381_final_exponent = bigint<12*bls12_381_q_limbs>("322277361516934140462891564586510139908379969514828494218366688025288661041104682794998680497580008899973249814104447692778988208376779573819485263026159588510513834876303014016798809919343532899164848730280942609956670917565618115867287399623286813270357901731510188149934363360381614501334086825442271920079363289954510565375378443704372994881406797882676971082200626541916413184642520269678897559532260949334760604962086348898118982248842634379637598665468817769075878555493752214492790122785850202957575200176084204422751485957336465472324810982833638490904279282696134323072515220044451592646885410572234451732790590013479358343841220074174848221722017083597872017638514103174122784843925578370430843522959600095676285723737049438346544753168912974976791528535276317256904336520179281145394686565050419250614107803233314658825463117900250701199181529205942363159325765991819433914303908860460720581408201373164047773794825411011922305820065611121544561808414055302212057471395719432072209245600258134364584636810093520285711072578721435517884103526483832733289802426157301542744476740008494780363354305116978805620671467071400711358839553375340724899735460480144599782014906586543813292157922220645089192130209334926661588737007768565838519456601560804957985667880395221049249803753582637708560");
    bls12_381_final_exponent_z = bigint<bls12_381_q_limbs>("15132376222941642752");
    bls12_381_final_exponent_is_z_neg = true;
}

} // namespace libff
