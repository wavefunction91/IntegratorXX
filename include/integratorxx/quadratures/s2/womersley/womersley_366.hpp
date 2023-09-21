#pragma once

namespace IntegratorXX {
namespace WomersleyGrids {

/**
 *  \brief Womersley Quadrature specification for index 26 grid with 366 points
 * 
 */
template <typename T>
struct womersley_366 {

  static constexpr std::array<cartesian_pt_t<T>,366> points = {
     0.0000000000000000e+00,      0.0000000000000000e+00,      1.0000000000000000e+00,
     1.8083162179353732e-01,      0.0000000000000000e+00,      9.8351406932464325e-01,
     1.1319022205603861e-01,      1.7152342998831921e-01,      9.7865606144137629e-01,
    -3.8861610343204810e-02,      2.5555569187134680e-01,      9.6601297278747256e-01,
    -1.5293014246514647e-01,      1.2098830710805794e-01,      9.8080283496160237e-01,
    -1.9171427325199314e-01,     -6.0962823631233051e-02,      9.7955559901741529e-01,
    -1.9279370766771053e-02,     -1.8188881270496313e-01,      9.8313008583585526e-01,
     1.6194793787063502e-01,     -2.0160885084420133e-01,      9.6598485323566519e-01,
     3.2297793988003776e-01,     -1.3583545744669528e-01,      9.3660769741183192e-01,
     3.5422848310436533e-01,      7.0031957361404759e-02,      9.3253295207499809e-01,
     3.0477077162315441e-01,      2.2702903983578410e-01,      9.2497167082860920e-01,
     1.7076608724925965e-01,      3.5152700771345236e-01,      9.2047145870668079e-01,
    -3.7859381168612346e-03,      4.3447528859691098e-01,      9.0067579642799678e-01,
    -1.8435958059207394e-01,      3.7736334387777387e-01,      9.0752876083422263e-01,
    -2.9529827323390229e-01,      2.2671955923205728e-01,      9.2811484811239675e-01,
    -3.5914685215033132e-01,      4.7934876616468038e-02,      9.3204924021978153e-01,
    -3.7026166463780891e-01,     -1.3714867307508424e-01,      9.1874726730117817e-01,
    -2.1796061921817639e-01,     -2.4766431988639817e-01,      9.4401035647138798e-01,
    -4.8522387040829597e-02,     -3.4991198764044307e-01,      9.3552508189859573e-01,
     1.5349420140316572e-01,     -3.7356443672939627e-01,      9.1481645249014476e-01,
     3.2667133630320383e-01,     -3.5401787864537548e-01,      8.7633166075254132e-01,
     4.5891116926430492e-01,     -2.4172274388990739e-01,      8.5496821801210998e-01,
     4.9113522841931806e-01,     -4.5261869782004763e-02,      8.6990663323677508e-01,
     5.2117331307519110e-01,      1.5407888164324682e-01,      8.3942723089604099e-01,
     4.7245466055127844e-01,      3.0934720492261147e-01,      8.2528231565322185e-01,
     3.3932169742533602e-01,      4.2465534313910952e-01,      8.3936203464285530e-01,
     1.6210964011517684e-01,      5.2542936629152337e-01,      8.3525112727862627e-01,
    -3.0976878634095779e-02,      5.9362444679010673e-01,      8.0414578850058538e-01,
    -2.0260935695089699e-01,      5.4559908732683782e-01,      8.1318576253157904e-01,
    -3.5840562569474532e-01,      4.0919845150164519e-01,      8.3910788028656591e-01,
    -4.6616581352287129e-01,      2.2945268788432790e-01,      8.5442430813104608e-01,
    -5.1987623832320551e-01,      3.1556807408428963e-02,      8.5365851763635614e-01,
    -5.4710131694651887e-01,     -1.5376400749890409e-01,      8.2282244682146444e-01,
    -4.4787536189581029e-01,     -3.0076890583259625e-01,      8.4198914808384540e-01,
    -2.9037117727799217e-01,     -3.9204023900512946e-01,      8.7291983045809762e-01,
    -1.3376144628300335e-01,     -4.9409163471068157e-01,      8.5905839847894261e-01,
     5.3096559630590630e-02,     -5.1771093870933271e-01,      8.5390639961069315e-01,
     2.4024166615536535e-01,     -5.4014409584693668e-01,      8.0655334452507810e-01,
     4.2283718608010468e-01,     -4.7319257860855490e-01,      7.7285024268460034e-01,
     5.7226923810430685e-01,     -3.7454970848461289e-01,      7.2953439603188530e-01,
     6.1727485517660319e-01,     -2.0639285338524102e-01,      7.5919282348965988e-01,
     6.3600554961083267e-01,     -1.4566020537579586e-02,      7.7154699915813396e-01,
     6.7037869516865567e-01,      1.7156875596265431e-01,      7.2191174463461794e-01,
     6.4626897101117686e-01,      3.3561596021406614e-01,      6.8534541973937868e-01,
     5.2949794949286011e-01,      4.6316247991575754e-01,      7.1071262735450402e-01,
     3.6436221717451522e-01,      5.6923433058421391e-01,      7.3702947809433861e-01,
     1.8987584205365263e-01,      6.5995817887759889e-01,      7.2691290175438461e-01,
     2.3795455164473265e-02,      7.2456972718191237e-01,      6.8879059718106250e-01,
    -1.7353913546652242e-01,      7.0902958088862000e-01,      6.8349193256865870e-01,
    -3.4198184947843069e-01,      6.0738360539052549e-01,      7.1703108058864395e-01,
    -4.8227955806420442e-01,      4.6121880661600057e-01,      7.4477086429123140e-01,
    -5.8440380904022649e-01,      2.7779926958743684e-01,      7.6243016322543344e-01,
    -6.5100468518595367e-01,      7.3587013473373777e-02,      7.5549841251587480e-01,
    -6.9190874328362395e-01,     -1.1529853367855172e-01,      7.1271911655241305e-01,
    -6.4176249387357287e-01,     -2.9962628328334601e-01,      7.0594970913159272e-01,
    -5.2386059734747237e-01,     -4.3483831653145005e-01,      7.3245185031020565e-01,
    -3.6803724254740977e-01,     -5.2416936656058732e-01,      7.6798116074391543e-01,
    -2.1783591161138469e-01,     -6.2442852777902513e-01,      7.5009101401640332e-01,
    -4.6220776803969389e-02,     -6.6578505748141426e-01,      7.4471061294042773e-01,
     1.2895833598571477e-01,     -6.7124370449468207e-01,      7.2993262480591414e-01,
     3.4781015089730805e-01,     -6.3847311609336288e-01,      6.8657132110132812e-01,
     5.4668244465582738e-01,     -5.4378975498967175e-01,      6.3673464415987369e-01,
     6.9975013342042180e-01,     -4.2620561506126381e-01,      5.7332235650491736e-01,
     7.4313069299671619e-01,     -2.6788809680603942e-01,      6.1318246934812015e-01,
     7.5794089251324959e-01,     -8.7639252604506460e-02,      6.4640928587013824e-01,
     7.8828989561666729e-01,      8.3887708452076173e-02,      6.0955876897910621e-01,
     7.9982074223861821e-01,      2.5306878192131893e-01,      5.4428207016373953e-01,
     7.3064420951848663e-01,      4.1934017346655655e-01,      5.3880688378503849e-01,
     5.8163206069743123e-01,      5.6698708972758527e-01,      5.8329219611709437e-01,
     4.1540783648056823e-01,      6.7453606347027761e-01,      6.1027651803797522e-01,
     2.5536547222466383e-01,      7.6807632684147098e-01,      5.8723694684615035e-01,
     1.0589102339592359e-01,      8.3428830031817358e-01,      5.4106388080925982e-01,
    -8.5840086910928745e-02,      8.2777635552916373e-01,      5.5445268933064051e-01,
    -2.8988781218846676e-01,      7.7450523533797810e-01,      5.6223366741831415e-01,
    -4.5449170208680573e-01,      6.5783983282138425e-01,      6.0056976870949097e-01,
    -5.9347069106236972e-01,      5.0798499582322587e-01,      6.2429462825530579e-01,
    -6.8764253279494614e-01,      3.4403741359041096e-01,      6.3936375025596492e-01,
    -7.4649964463934348e-01,      1.5363426243270881e-01,      6.4740620475872135e-01,
    -8.0132782397623337e-01,     -4.7118154826680692e-02,      5.9636699942840887e-01,
    -7.8569881577064626e-01,     -2.4514459840782726e-01,      5.6796258395079968e-01,
    -7.0048489083429177e-01,     -4.2000414061834696e-01,      5.7698998221486841e-01,
    -5.6135152243854547e-01,     -5.6309901448412947e-01,      6.0646843952750684e-01,
    -4.0623494231340418e-01,     -6.5822346580703783e-01,      6.3380994052207507e-01,
    -2.7046256752916137e-01,     -7.5583770550655627e-01,      5.9628798621145995e-01,
    -1.0286838735552602e-01,     -7.8777779627598588e-01,      6.0730901407555660e-01,
     7.6276377339914270e-02,     -8.0544019305926673e-01,      5.8774825364653382e-01,
     2.4857567876170983e-01,     -7.6591970270257870e-01,      5.9293940747782048e-01,
     4.5669623085815353e-01,     -6.9472646299161911e-01,      5.5568308804489508e-01,
     6.3254060624775998e-01,     -6.0510259486687545e-01,      4.8346998989915624e-01,
     7.7326631208404528e-01,     -4.8155429784592746e-01,      4.1251020450657422e-01,
     8.4653079983104151e-01,     -3.0500516306460140e-01,      4.3629973119560023e-01,
     8.6342546938303166e-01,     -1.3819388370189015e-01,      4.8517925484101215e-01,
     8.9063727245004665e-01,      2.8424794940502447e-02,      4.5382516452410582e-01,
     8.9572177415863119e-01,      2.0013291138782421e-01,      3.9702559247175429e-01,
     8.4272052924739804e-01,      3.8802959949534688e-01,      3.7316904949428220e-01,
     7.0840577417699280e-01,      5.5838096846141405e-01,      4.3170818057200261e-01,
     5.5090686989099991e-01,      6.9949731901030221e-01,      4.5519789257453780e-01,
     3.9282128225628760e-01,      8.0324867921329313e-01,      4.4775327977428153e-01,
     2.5775295129637915e-01,      8.8693518609394595e-01,      3.8328760972212983e-01,
     7.9337243601899129e-02,      9.2145818808347235e-01,      3.8028990177439220e-01,
    -1.0988719682575672e-01,      9.0903492230658556e-01,      4.0197053872247573e-01,
    -2.9141523674734032e-01,      8.6365059990109361e-01,      4.1132080069207849e-01,
    -4.6414626265100595e-01,      7.6585990466846587e-01,      4.4500208234154576e-01,
    -6.1038483016164946e-01,      6.3194776732681868e-01,      4.7756923946081742e-01,
    -7.3523772146377497e-01,      4.7836860459728098e-01,      4.8019680452123681e-01,
    -8.1015267223559706e-01,      3.0173260733542190e-01,      5.0260330414760435e-01,
    -8.5769766313224394e-01,      1.1079838522164953e-01,      5.0207413445602112e-01,
    -8.9528894676720450e-01,     -8.4225854248724780e-02,      4.3745137703811393e-01,
    -8.6720531928121514e-01,     -2.7027202770779934e-01,      4.1821999623294004e-01,
    -7.8693994675971368e-01,     -4.4998729528816978e-01,      4.2218118654560682e-01,
    -6.5888302751040773e-01,     -5.9688435312413490e-01,      4.5782335573264743e-01,
    -5.0398341755283216e-01,     -7.2335483880636497e-01,      4.7197297804766891e-01,
    -3.7666802085140683e-01,     -8.2096108965779646e-01,      4.2912013625063961e-01,
    -2.0495304111818008e-01,     -8.7656519751072481e-01,      4.3546263381534078e-01,
    -2.5762118949995121e-02,     -8.9713742576402289e-01,      4.4099971941114513e-01,
     1.5224643764093898e-01,     -8.9821283077434066e-01,      4.1235268018771160e-01,
     3.1937473185884308e-01,     -8.3854668691135992e-01,      4.4140597472176774e-01,
     4.8410026213346297e-01,     -7.7832819351352556e-01,      3.9981515652146588e-01,
     6.4073221765609689e-01,     -6.9061990819275831e-01,      3.3541968884566981e-01,
     7.7531918447795667e-01,     -5.7564727982403896e-01,      2.5982758015965274e-01,
     8.8169563710793719e-01,     -3.9257296784797957e-01,      2.6172364895029604e-01,
     9.3548567657774528e-01,     -2.0625791795148088e-01,      2.8692197580561590e-01,
     9.5716418180246787e-01,     -3.5463173995822612e-02,      2.8736578147816111e-01,
     9.6377114160522781e-01,      1.3650439607574757e-01,      2.2915439437407975e-01,
     9.2393912957966373e-01,      3.1620758355492778e-01,      2.1528875710061326e-01,
     8.2438946366103882e-01,      5.0737984982233975e-01,      2.5089380263156136e-01,
     6.9812800114496987e-01,      6.5656510487947417e-01,      2.8555132125755256e-01,
     5.5111060541709778e-01,      7.8104713881559107e-01,      2.9367067872836522e-01,
     4.0470195283566274e-01,      8.8106190402504436e-01,      2.4483923428806195e-01,
     2.2857693314802199e-01,      9.5350602127604112e-01,      1.9641500203135837e-01,
     5.0157449022563978e-02,      9.7733834778438677e-01,      2.0565501272187325e-01,
    -1.3075939895087774e-01,      9.6451153982973936e-01,      2.2938933959813865e-01,
    -3.1084444453129229e-01,      9.2145542065323716e-01,      2.3301424645887681e-01,
    -4.6640408107614512e-01,      8.4258409950215341e-01,      2.6929401854044555e-01,
    -6.1175454695450648e-01,      7.2948710318303589e-01,      3.0594924508831417e-01,
    -7.4133058762250692e-01,      5.8774017001685108e-01,      3.2402230232474960e-01,
    -8.4887205529464849e-01,      4.0932886163546367e-01,      3.3446392447027046e-01,
    -9.1032758371347300e-01,      2.2103064923851612e-01,      3.4992733878276761e-01,
    -9.4822866363307434e-01,      4.3451915089680888e-02,      3.1460186353496589e-01,
    -9.5870749076197226e-01,     -1.3503921858606241e-01,      2.5028854668272099e-01,
    -9.1824180891342444e-01,     -3.1041772070781515e-01,      2.4591221814697856e-01,
    -8.3893724120148128e-01,     -4.7913868740366178e-01,      2.5808995245523147e-01,
    -7.2018063546453781e-01,     -6.3136693559336365e-01,      2.8760327700034999e-01,
    -5.8607580896525069e-01,     -7.5093328153177297e-01,      3.0432606334924733e-01,
    -4.5649812876262552e-01,     -8.5668105481263834e-01,      2.4022287310188983e-01,
    -2.8304748338977587e-01,     -9.2227171766288152e-01,      2.6324703406852029e-01,
    -9.0284124201099630e-02,     -9.5918965607302387e-01,      2.6796264776971152e-01,
     9.2515096107936998e-02,     -9.6559451183828537e-01,      2.4303949411550896e-01,
     2.8679519686487259e-01,     -9.2309240655182034e-01,      2.5622045980289687e-01,
     4.4916543235793976e-01,     -8.6252564348084948e-01,      2.3302345099292360e-01,
     5.9795048720831956e-01,     -7.8135453809699651e-01,      1.7871849552456046e-01,
     7.2901846608859877e-01,     -6.7553912441420216e-01,      1.1035835939121089e-01,
     8.5771168138853027e-01,     -5.0449749629578100e-01,      9.9060324252188470e-02,
     9.4328425840989583e-01,     -3.1599526852829091e-01,      1.0179291774886082e-01,
     9.8610846926891116e-01,     -1.2703179359376071e-01,      1.0701873785686862e-01,
     9.9578409920247835e-01,      5.9662398931442991e-02,      6.9673710459928873e-02,
     9.6764475525249549e-01,      2.5095745336975001e-01,      2.6153092178704696e-02,
     9.0891234019330580e-01,      4.1355721865162859e-01,      5.3374008144951871e-02,
     8.1715636046296469e-01,      5.6988430929491019e-01,      8.6529512736316372e-02,
     6.8834599727803847e-01,      7.1561825834306236e-01,      1.1861828003029216e-01,
     5.3838593748685248e-01,      8.3619514086077151e-01,      1.0449051974814462e-01,
     3.5853171673112605e-01,      9.3183411618659517e-01,      5.6037380457856910e-02,
     1.6924558267295114e-01,      9.8555250644483483e-01,      6.4953665021492430e-03,
    -4.0073270058352418e-03,      9.9978665538633094e-01,      2.0262898156045907e-02,
    -1.8160332471212931e-01,      9.8214946154276483e-01,      4.9017013829465661e-02,
    -3.5835053545211504e-01,      9.3248772767311616e-01,      4.5293834902895944e-02,
    -5.0948676567983520e-01,      8.5694978343564920e-01,      7.7847956085530329e-02,
    -6.4421267345046218e-01,      7.5446163817372258e-01,      1.2560918712431524e-01,
    -7.7179077419378539e-01,      6.1699988208695444e-01,      1.5378604089462003e-01,
    -8.7459207696098928e-01,      4.5350227681976246e-01,      1.7153537196845015e-01,
    -9.4645741495951785e-01,      2.7036101402962281e-01,      1.7641792358210368e-01,
    -9.8568436351899913e-01,      1.0033570734881309e-01,      1.3549568755115007e-01,
    -9.9477075446034691e-01,     -6.3446956136449364e-02,      8.0035178686698877e-02,
    -9.6577913254000691e-01,     -2.5218662706736106e-01,      6.0601751448777580e-02,
    -9.0181818957893267e-01,     -4.2683967678429086e-01,      6.7318966623518967e-02,
    -8.0606907120053528e-01,     -5.8429641048008485e-01,      9.4076336844046385e-02,
    -6.8646116455543149e-01,     -7.1825966172834055e-01,      1.1346421414300881e-01,
    -5.4827583243373135e-01,     -8.3153928602051774e-01,      8.9084383443936613e-02,
    -3.4043631100207594e-01,     -9.3628655723198284e-01,      8.6432637920976454e-02,
    -1.3995347266436731e-01,     -9.8482932397193279e-01,      1.0258766073056735e-01,
     7.4827121665631019e-02,     -9.9433328006415223e-01,      7.5513111577403008e-02,
     2.7227740657782085e-01,     -9.5873238999955557e-01,      8.1836533608140416e-02,
     4.4922956168445888e-01,     -8.9262742885937851e-01,      3.7537636535672145e-02,
     5.8213065947708142e-01,     -8.1279689261015708e-01,     -2.2025137004132802e-02,
     7.2120587116454116e-01,     -6.9001331832768420e-01,     -6.1185880137603121e-02,
     8.4588374554161461e-01,     -5.3004850343413690e-01,     -5.9407684988735360e-02,
     9.4063871603793148e-01,     -3.3351201389801222e-01,     -6.2996368754110674e-02,
     9.9005761469842679e-01,     -1.2500080342114428e-01,     -6.4503633399444746e-02,
     9.9024113432303451e-01,      9.0716036082702997e-02,     -1.0579743235103430e-01,
     9.4744338078627399e-01,      2.6984034277774344e-01,     -1.7186398579651693e-01,
     8.9336225010297576e-01,      4.2434260648251043e-01,     -1.4777429551372392e-01,
     8.1117464572264719e-01,      5.7720538781807562e-01,     -9.3966134381078106e-02,
     6.8691877981244343e-01,      7.2473232181831415e-01,     -5.3904096809235066e-02,
     5.2183983625377861e-01,      8.5051746028929387e-01,     -6.5599047566858026e-02,
     3.3757314323703541e-01,      9.3356938157858416e-01,     -1.2038514336930499e-01,
     1.5972066301316770e-01,      9.6976780202902524e-01,     -1.8449856355654257e-01,
    -1.6487401809991902e-02,      9.8445050251072030e-01,     -1.7488674531806636e-01,
    -1.9079520053036905e-01,      9.7176650722574820e-01,     -1.3880578838379223e-01,
    -3.6203972270807605e-01,      9.2048725567249456e-01,     -1.4707294559496309e-01,
    -5.1820189197302635e-01,      8.4594619060322473e-01,     -1.2586437843754070e-01,
    -6.5704104408467079e-01,      7.5082820850405474e-01,     -6.7483832898831433e-02,
    -7.8153616772134127e-01,      6.2326987446312154e-01,     -2.7127147475624769e-02,
    -8.8398291829076447e-01,      4.6740832419887945e-01,     -1.0181288707177161e-02,
    -9.5374810910285479e-01,      3.0053337343254138e-01,     -6.6510026301159052e-03,
    -9.8848295141192966e-01,      1.3757445212585681e-01,     -6.3045419264456226e-02,
    -9.9356860805347769e-01,     -3.1854440771006416e-02,     -1.0865871200157534e-01,
    -9.6978454196331354e-01,     -2.1781789125954748e-01,     -1.0987860763701987e-01,
    -9.0059980696729491e-01,     -4.2092975521159814e-01,     -1.0834264565707781e-01,
    -8.0125815370483688e-01,     -5.9091035322228358e-01,     -9.3863334568041287e-02,
    -6.8258420173848267e-01,     -7.2662040500503133e-01,     -7.8112704263536337e-02,
    -5.2221449082570048e-01,     -8.4913999670286822e-01,     -7.9077756487569678e-02,
    -3.5470613151350189e-01,     -9.3264170016495385e-01,     -6.6053155716821396e-02,
    -1.2911275419512935e-01,     -9.9011660282107239e-01,     -5.4763213220253161e-02,
     7.9553486124642692e-02,     -9.9305563720971202e-01,     -8.6670319322297115e-02,
     2.8728761408466336e-01,     -9.5347418242763038e-01,     -9.1393710054375008e-02,
     4.3692650067458755e-01,     -8.8211608068737579e-01,     -1.7597287632189995e-01,
     5.7792811167424218e-01,     -7.8582243362272086e-01,     -2.2018674018184931e-01,
     7.1347826368724354e-01,     -6.6019167796697731e-01,     -2.3472476560640809e-01,
     8.3615529741050976e-01,     -5.0142159041048451e-01,     -2.2230768606272364e-01,
     9.2241569675702439e-01,     -3.1840536658304786e-01,     -2.1855732636397290e-01,
     9.6895656861319235e-01,     -1.0917329353371517e-01,     -2.2182055838075351e-01,
     9.5671677167709035e-01,      7.7547014406234321e-02,     -2.8049862628620015e-01,
     9.0078655874990021e-01,      2.3537271718214478e-01,     -3.6494281686561048e-01,
     8.4822190905784012e-01,      4.0325051025583858e-01,     -3.4337824475740908e-01,
     7.7536395932903912e-01,      5.6844058106883155e-01,     -2.7511095283126275e-01,
     6.6095577510426440e-01,      7.1679298725782026e-01,     -2.2213796788106979e-01,
     5.0200059963737376e-01,      8.3260510023430712e-01,     -2.3401740325782763e-01,
     3.3048206011601572e-01,      8.9377049153332400e-01,     -3.0324233940159889e-01,
     1.5803101788581225e-01,      9.1571015492651175e-01,     -3.6946056562296076e-01,
    -2.4199958872727134e-02,      9.3364730558613285e-01,     -3.5737525202833043e-01,
    -2.0161876747928553e-01,      9.2284478471176190e-01,     -3.2818801917564788e-01,
    -3.6308290460492293e-01,      8.6842134263772863e-01,     -3.3766133334294751e-01,
    -5.3016322529462756e-01,      7.9217020769788260e-01,     -3.0231327556210896e-01,
    -6.8375095886177328e-01,      6.8856947288031356e-01,     -2.4157132957562033e-01,
    -8.0593723760676084e-01,      5.5577542759259835e-01,     -2.0390890888592389e-01,
    -8.9690958967680123e-01,      3.9749670666983489e-01,     -1.9377707844950912e-01,
    -9.4598956549285651e-01,      2.3283370161296604e-01,     -2.2559301711675522e-01,
    -9.5650595647786907e-01,      6.0695666919224144e-02,     -2.8532856716352029e-01,
    -9.4922134896320665e-01,     -1.2132727297547195e-01,     -2.9027318771255051e-01,
    -9.1084288624366261e-01,     -3.2032219398147710e-01,     -2.6030545254029419e-01,
    -8.0936183519125138e-01,     -5.2660448207332478e-01,     -2.6004064912266123e-01,
    -6.6134247735158835e-01,     -7.0780268803741653e-01,     -2.4827702764748769e-01,
    -5.1179468136945316e-01,     -8.1987553681281899e-01,     -2.5665211524129616e-01,
    -3.5492224034414638e-01,     -8.9516194110489100e-01,     -2.6965775068856351e-01,
    -2.0482156952850900e-01,     -9.5692718625717210e-01,     -2.0576317663228669e-01,
    -6.3456860963803808e-03,     -9.6993771272114537e-01,     -2.4327056071222228e-01,
     1.9933733580620125e-01,     -9.4754141255765512e-01,     -2.4985975674751909e-01,
     3.5061066850853412e-01,     -8.7359540893014509e-01,     -3.3749551200596967e-01,
     5.0027711401595898e-01,     -7.7101995745139662e-01,     -3.9401907873034470e-01,
     6.3898738260921084e-01,     -6.5181213953971606e-01,     -4.0845570092102723e-01,
     7.6576647880319371e-01,     -5.1171597333582419e-01,     -3.8954905028035669e-01,
     8.5188832468687237e-01,     -3.4723430604918570e-01,     -3.9206455968977133e-01,
     9.0814192347283806e-01,     -1.7589860100870922e-01,     -3.7991831884529198e-01,
     8.9957151433652571e-01,     -1.2268139371886706e-03,     -4.3677177738706230e-01,
     8.4133744517438347e-01,      1.5950275793310545e-01,     -5.1643990314379684e-01,
     7.7738288348376428e-01,      3.7097750072563446e-01,     -5.0798774239328859e-01,
     7.0979507956446575e-01,      5.5334329637375534e-01,     -4.3589235068340004e-01,
     5.9411666321837187e-01,      7.0345466506893073e-01,     -3.9009860891706083e-01,
     4.4073611335469026e-01,      7.9102405471426818e-01,     -4.2430251383700307e-01,
     2.7084289636148695e-01,      8.2073762949534423e-01,     -5.0302472008926546e-01,
     7.3068908047543554e-02,      8.4388527076677433e-01,     -5.3152477313821178e-01,
    -1.1593340192078976e-01,      8.5075206492244948e-01,     -5.1262107872117479e-01,
    -2.7930118334928572e-01,      8.0994597848172045e-01,     -5.1573089971512964e-01,
    -4.5829186660218385e-01,      7.4990449465143627e-01,     -4.7708679913393132e-01,
    -6.3363193846039112e-01,      6.5357699401570390e-01,     -4.1394163774175174e-01,
    -7.6917022713950289e-01,      5.1369025576777028e-01,     -3.8013087589856259e-01,
    -8.5410742865947187e-01,      3.5437353500337831e-01,     -3.8068346168163730e-01,
    -8.8720413519584385e-01,      1.7805771256904154e-01,     -4.2563396655586044e-01,
    -8.8244698091658469e-01,     -2.6815047977068970e-03,     -4.7040422553716998e-01,
    -8.7814476657403495e-01,     -1.9234558117079600e-01,     -4.3802391069746666e-01,
    -8.2891394155208820e-01,     -3.8968059018212370e-01,     -4.0131124471648305e-01,
    -7.0839901669356653e-01,     -5.9126712778435875e-01,     -3.8545300199793858e-01,
    -5.5560721749259367e-01,     -7.1130150487679278e-01,     -4.3052385419410577e-01,
    -4.0996267046818174e-01,     -7.9271963393933831e-01,     -4.5113877110006678e-01,
    -2.3991679451084605e-01,     -8.7170533135664374e-01,     -4.2728181215217254e-01,
    -8.4268380449495223e-02,     -9.1743406366604441e-01,     -3.8886190207016641e-01,
     1.2329714755030069e-01,     -9.0429112395282862e-01,     -4.0872408363844820e-01,
     2.8026721774044977e-01,     -8.2283355492247146e-01,     -4.9436345693596268e-01,
     4.3827411244727493e-01,     -7.1243844458131589e-01,     -5.4803947398076058e-01,
     5.9389064699281513e-01,     -5.7378054144851121e-01,     -5.6397676341273884e-01,
     7.0849100359408868e-01,     -4.3460218062299227e-01,     -5.5602287940693684e-01,
     7.8302111518862350e-01,     -2.7413401292205447e-01,     -5.5832649599317352e-01,
     8.1728217354701438e-01,     -1.0922129327638352e-01,     -5.6579197404814952e-01,
     7.6367359574550664e-01,      5.6006727780332125e-02,     -6.4316862921357443e-01,
     7.2415662235663258e-01,      2.6211014332860999e-01,     -6.3788357798370143e-01,
     6.4415979982676530e-01,      4.6603638459927993e-01,     -6.0652142626355221e-01,
     5.5412170515037851e-01,      6.1849896761861867e-01,     -5.5714285684727205e-01,
     3.9844566728517844e-01,      6.9704893011454439e-01,     -5.9612401331253051e-01,
     1.9939987907566312e-01,      7.3204175764855128e-01,     -6.5142501739143455e-01,
    -1.0087511613910602e-02,      7.4608848946334039e-01,     -6.6577038684500700e-01,
    -1.9055835450796929e-01,      7.1435424454963747e-01,     -6.7333908754886040e-01,
    -3.7386576076760397e-01,      6.7360166607615934e-01,     -6.3756190945278646e-01,
    -5.4909604451143668e-01,      6.0523054344803562e-01,     -5.7635884931133818e-01,
    -6.8912492025074434e-01,      4.8298312925618786e-01,     -5.4021675385284540e-01,
    -7.7558473966059160e-01,      3.1221068191066442e-01,     -5.4862810874625290e-01,
    -7.8786772573914376e-01,      1.3546851470221266e-01,     -6.0076012539365964e-01,
    -7.8778141300413740e-01,     -6.2968499073392178e-02,     -6.1272784615165710e-01,
    -7.7915900463263688e-01,     -2.5967113542034870e-01,     -5.7051042666141116e-01,
    -7.2105745757597439e-01,     -4.4740227021683998e-01,     -5.2906271034622121e-01,
    -5.9056961131312602e-01,     -5.8002486356648975e-01,     -5.6106923978965229e-01,
    -4.3745216717069496e-01,     -6.4769489775015499e-01,     -6.2380038543277494e-01,
    -2.9188149770269384e-01,     -7.4848105377479146e-01,     -5.9546729837918888e-01,
    -1.1221687118047291e-01,     -8.0882666738345377e-01,     -5.7724067420084058e-01,
     5.6466302754157809e-02,     -8.3440397802941535e-01,     -5.4825318795421751e-01,
     2.4386733194910959e-01,     -7.3139780511653618e-01,     -6.3685632216280630e-01,
     4.3408744224301438e-01,     -5.9311952341938401e-01,     -6.7806881909263494e-01,
     5.6521401428567752e-01,     -4.2910689223915444e-01,     -7.0455687711349790e-01,
     6.4858639090262460e-01,     -2.7463731605807584e-01,     -7.0986621145418971e-01,
     6.8310319030380218e-01,     -1.1024023378043590e-01,     -7.2195368428508089e-01,
     6.3181574800971674e-01,      7.8358516506991785e-02,     -7.7114771831196238e-01,
     5.9882110818732937e-01,      2.5652428819006429e-01,     -7.5868871743151656e-01,
     5.3221703791545194e-01,      4.2661257148217985e-01,     -7.3126379536106145e-01,
     4.1439700251748857e-01,      5.6216096206091448e-01,     -7.1571654797082040e-01,
     2.3443903241769731e-01,      6.0755723521080662e-01,     -7.5888902088648935e-01,
     3.8408778666278301e-02,      6.2390582514669557e-01,     -7.8055511469042715e-01,
    -1.7257274880229825e-01,      5.8869552305467954e-01,     -7.8971908138666380e-01,
    -3.6487971652936524e-01,      5.4186599110141764e-01,     -7.5712881344796834e-01,
    -5.2273074761171490e-01,      4.7503873206448832e-01,     -7.0787765082665288e-01,
    -6.4200081873636516e-01,      3.4974101580030847e-01,     -6.8228745452983752e-01,
    -6.6947369150794156e-01,      1.6394588519837364e-01,     -7.2451826968355404e-01,
    -6.7172749640005402e-01,     -3.5905505144554685e-02,     -7.3992767570920726e-01,
    -6.7027985920758193e-01,     -2.3723887525129259e-01,     -7.0316614424342549e-01,
    -6.1024329293964430e-01,     -4.1085246409387077e-01,     -6.7735026106894958e-01,
    -4.5871007579823464e-01,     -4.9753329030457027e-01,     -7.3623752376518148e-01,
    -2.8120992121451238e-01,     -6.0082042066115726e-01,     -7.4828858225091055e-01,
    -1.2460775698856802e-01,     -6.8775247067333933e-01,     -7.1517092081613287e-01,
     7.9333174793371680e-02,     -7.1400702601394228e-01,     -6.9562936552443633e-01,
     2.7030846932501290e-01,     -5.9248842625423270e-01,     -7.5887469068743574e-01,
     4.1137803951797319e-01,     -4.3486093209813864e-01,     -8.0103937377452739e-01,
     4.9123275269605698e-01,     -2.5505431724934602e-01,     -8.3284913275522143e-01,
     5.2887962655851528e-01,     -8.5735588454678893e-02,     -8.4435522707190924e-01,
     4.6978765215850199e-01,      8.7617674769826598e-02,     -8.7842057406878349e-01,
     4.3538713783119726e-01,      2.6802941076521020e-01,     -8.5941740451076054e-01,
     3.4753713850901219e-01,      4.2264815339317119e-01,     -8.3701043947510601e-01,
     1.8090577333044139e-01,      4.6058671586139505e-01,     -8.6898387692046941e-01,
    -1.1284676517849420e-02,      4.8685950523895050e-01,     -8.7340739533986766e-01,
    -2.1887764232486040e-01,      4.3541073925627771e-01,     -8.7321822348746947e-01,
    -3.9014940914367063e-01,      3.6338849170875454e-01,     -8.4600959961366928e-01,
    -5.1568530112410249e-01,      2.6627753597662490e-01,     -8.1434939923767424e-01,
    -5.3801577862390859e-01,      7.5709904080929077e-02,     -8.3952786277512303e-01,
    -5.4039989842742719e-01,     -1.3199139834912374e-01,     -8.3099110737809301e-01,
    -5.0108246742932161e-01,     -3.0512428805326541e-01,     -8.0982438199583173e-01,
    -3.3403792451110714e-01,     -4.1770879473085432e-01,     -8.4494853558829774e-01,
    -1.5770121106503218e-01,     -5.0711373234679125e-01,     -8.4732873814944409e-01,
     2.0341907418782435e-02,     -5.9201568058855314e-01,     -8.0566968463498601e-01,
     1.7407581749403106e-01,     -4.9679977892101979e-01,     -8.5022796321210858e-01,
     2.9432302548571421e-01,     -3.4544235903381465e-01,     -8.9109120366777761e-01,
     3.4946402882671407e-01,     -1.6151361595232341e-01,     -9.2292374788939469e-01,
     3.2405027602296216e-01,      2.1125675248532684e-02,     -9.4580395667111461e-01,
     2.7553632434534026e-01,      2.2913300160191263e-01,     -9.3358331258820038e-01,
     1.2152375981103457e-01,      3.0609028962489543e-01,     -9.4421433498900964e-01,
    -6.8601007497540750e-02,      3.2290840188983483e-01,     -9.4394071093436549e-01,
    -2.4360573134797164e-01,      2.4574206708580568e-01,     -9.3822549747851935e-01,
    -3.7215887107441098e-01,      1.3818597102293201e-01,     -9.1782482647456720e-01,
    -3.9114762795926333e-01,     -4.3118741320380231e-02,     -9.1931730500887876e-01,
    -3.5113311153077370e-01,     -2.2712163583422462e-01,     -9.0836187751507003e-01,
    -1.8763957246692706e-01,     -3.1588207617272618e-01,     -9.3005908672365378e-01,
    -1.0862087269264378e-02,     -4.1621883498615286e-01,     -9.0919959110358417e-01,
     1.2298494577441864e-01,     -3.0555596830748277e-01,     -9.4419820659887921e-01,
     1.7687791803844854e-01,     -1.3417719160011246e-01,     -9.7504393919694232e-01,
     1.5709929901731706e-01,      8.5447368650816022e-02,     -9.8387934089446105e-01,
    -1.8480798355532024e-02,      1.5109727548358143e-01,     -9.8834613037820962e-01,
    -1.7622625934941116e-01,      6.6326544598748649e-02,     -9.8211256737570785e-01,
    -2.1153510272618370e-01,     -9.3871368469085600e-02,     -9.7285202702999185e-01,
    -5.3177783259298859e-02,     -2.0886282619779997e-01,     -9.7649805079185625e-01,
     1.5152217608325816e-02,     -4.8303853877100345e-02,     -9.9871775192101686e-01
};


static constexpr auto weights = detail::create_array<366, T>(4.0 * M_PI / 366.0);
};
}  // namespace WomersleyGrids
}  // namespace IntegratorXX