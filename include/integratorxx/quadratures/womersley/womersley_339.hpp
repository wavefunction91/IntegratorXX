#pragma once

namespace IntegratorXX {
namespace WomersleyGrids {

/**
 *  \brief Womersley Quadrature specification for index 25 grid with 339 points
 * 
 */
template <typename T>
struct womersley_339 {

  static constexpr std::array<cartesian_pt_t<T>,339> points = {
     0.0000000000000000e+00,      0.0000000000000000e+00,      1.0000000000000000e+00,
     1.9266881650578010e-01,      0.0000000000000000e+00,      9.8126384176034021e-01,
     1.5601748654593295e-01,      1.7771630769471264e-01,      9.7163545523578332e-01,
     1.0195311655880239e-03,      2.6038632265292222e-01,      9.6550397385587738e-01,
    -1.3005139104052102e-01,      1.2284432918414932e-01,      9.8386782977985510e-01,
    -2.1193189721729566e-01,     -6.8614005226295030e-02,      9.7487280669258780e-01,
    -6.6460118656043635e-02,     -1.9782329542160873e-01,      9.7798210434381649e-01,
     1.1405905558783562e-01,     -1.9548494508824227e-01,      9.7405141962950703e-01,
     3.0174460330706199e-01,     -1.8550770504612057e-01,      9.3516687588022773e-01,
     3.8460979213489238e-01,     -5.2200861307098263e-03,      9.2306449313942462e-01,
     3.5325361644342584e-01,      1.8566711838721395e-01,      9.1691853706827697e-01,
     2.7142662924834865e-01,      3.5362798602553114e-01,      8.9513956031135511e-01,
     1.0635850767962805e-01,      4.2007053794962024e-01,      9.0123726675658455e-01,
    -1.0229945982939084e-01,      4.2825535679852039e-01,      8.9784863417615490e-01,
    -2.2748320983321493e-01,      3.1094613766535201e-01,      9.2280219370945227e-01,
    -3.2113555821923223e-01,      1.3976616522192314e-01,      9.3666289149639115e-01,
    -4.1231976430508938e-01,     -4.1721323473258611e-03,      9.1102964017372778e-01,
    -3.2495098427798086e-01,     -1.9775273048654948e-01,      9.2482469441612858e-01,
    -1.9043055512287560e-01,     -3.4534902298979780e-01,      9.1895062761586144e-01,
    -2.2820368775633312e-02,     -4.0836029804814872e-01,      9.1253553232023310e-01,
     1.5513759070524136e-01,     -3.7844707552094625e-01,      9.1252952773037221e-01,
     3.3497986334516061e-01,     -3.7864959724263930e-01,      8.6279370284051371e-01,
     4.7110001129569978e-01,     -2.2452213724350625e-01,      8.5302672246817668e-01,
     5.4355801270875559e-01,     -2.6029044624999914e-02,      8.3896792290052313e-01,
     5.2127479815729993e-01,      1.8343905206541361e-01,      8.3344027919426156e-01,
     4.5520556415181207e-01,      3.7564084570380740e-01,      8.0726813971824674e-01,
     3.2843360805919064e-01,      5.2469981393008325e-01,      7.8537982553600005e-01,
     1.3939279620326883e-01,      5.7826669524959107e-01,      8.0385152704448504e-01,
    -4.7817866376437958e-02,      5.9061334388742870e-01,      8.0553667183891453e-01,
    -2.3236094713843258e-01,      5.7731716531468802e-01,      7.8276004041975966e-01,
    -3.4691556608775026e-01,      4.5046313061917392e-01,      8.2263756172374525e-01,
    -4.3859333990607252e-01,      2.8947253073246887e-01,      8.5078877292861366e-01,
    -5.4999572511717087e-01,      1.6217891945961585e-01,      8.1926961400734788e-01,
    -5.7294599797873647e-01,     -2.9887383103790615e-02,      8.1904800087135121e-01,
    -5.0178444300783509e-01,     -2.1499427460814738e-01,      8.3784833630021205e-01,
    -3.9489057102640357e-01,     -3.7912491243798491e-01,      8.3685466939208220e-01,
    -2.4324362401620639e-01,     -5.0774427254676835e-01,      8.2645525775532180e-01,
    -1.0923870244919212e-01,     -6.1068744960892285e-01,      7.8430079993416180e-01,
     5.1023513173340404e-02,     -5.7052149561725884e-01,      8.1969617794771743e-01,
     2.3686695833462273e-01,     -5.5170721155143454e-01,      7.9969569010433261e-01,
     4.2408279105864705e-01,     -5.2105397713347423e-01,      7.4071353385860095e-01,
     5.4079762079060389e-01,     -3.7799200694278079e-01,      7.5143860430150311e-01,
     6.3801188233646033e-01,     -1.9036672814014224e-01,      7.4612421674591412e-01,
     6.8764973784560079e-01,      1.5772251106870538e-02,      7.2587125176294065e-01,
     6.5839698181508688e-01,      2.0877209164352853e-01,      7.2313735077616492e-01,
     5.9603443631566899e-01,      3.9947734590280287e-01,      6.9653485256397263e-01,
     4.7480190956259427e-01,      5.7115129568696266e-01,      6.6958893666996744e-01,
     2.7923974321303180e-01,      6.8394108613293036e-01,      6.7398053125415092e-01,
     9.4211777649343825e-02,      7.3637037258539162e-01,      6.6998717549711495e-01,
    -7.9137766545696334e-02,      7.4119872470331216e-01,      6.6660457874540779e-01,
    -2.6203392985943502e-01,      7.2048405248221570e-01,      6.4205992689251723e-01,
    -4.1920227362698759e-01,      6.1858712855294384e-01,      6.6454451933229264e-01,
    -5.1757016632630026e-01,      4.6781412415675194e-01,      7.1642938812448009e-01,
    -6.2535627006297057e-01,      3.2465394178468171e-01,      7.0959802252866899e-01,
    -7.0913324130602795e-01,      1.4216017101648012e-01,      6.9059433233365874e-01,
    -7.2061869030518499e-01,     -3.9843109576660152e-02,      6.9218583473089224e-01,
    -6.4503847906024081e-01,     -2.3525053244802327e-01,      7.2703682679392023e-01,
    -5.3451674552886841e-01,     -4.3363645114492028e-01,      7.2543178658483543e-01,
    -3.9793472201843694e-01,     -5.7753403759093236e-01,      7.1281301365507133e-01,
    -2.7764812604874517e-01,     -6.9513050464940640e-01,      6.6310263127775460e-01,
    -1.1504521810675294e-01,     -7.6921414849483649e-01,      6.2854927535248439e-01,
     4.8081576089708922e-02,     -7.4544197395540068e-01,      6.6483413383054057e-01,
     2.1240433115191357e-01,     -6.9758034526842549e-01,      6.8429968727385271e-01,
     4.0106743862353550e-01,     -6.7077563565997422e-01,      6.2386293068342669e-01,
     5.7288754797600272e-01,     -5.7082709449606606e-01,      5.8818048723433736e-01,
     6.6804519706605203e-01,     -4.1614454409390911e-01,      6.1687870209454654e-01,
     7.5091520827292435e-01,     -2.3690834725145882e-01,      6.1644203700511291e-01,
     8.0299096893344046e-01,     -4.6114839377530734e-02,      5.9420444747621914e-01,
     8.0331761443986927e-01,      1.5206011063059341e-01,      5.7581119569321437e-01,
     7.6212629431667767e-01,      3.1805762580352631e-01,      5.6391742141855550e-01,
     6.7544825551965137e-01,      4.8866138588973773e-01,      5.5224967546914905e-01,
     5.4727614224453880e-01,      6.4713315046998265e-01,      5.3076125489030368e-01,
     3.8136372736880875e-01,      7.5347506276769893e-01,      5.3557169196530652e-01,
     1.9992806593663184e-01,      8.3362843025491717e-01,      5.1487125451083460e-01,
     2.7512631775785661e-02,      8.6923571535163435e-01,      4.9363177192103702e-01,
    -1.4996596670223208e-01,      8.5367898081262450e-01,      4.9874082101807504e-01,
    -3.3646982628593725e-01,      8.1159319736932167e-01,      4.7760290826475527e-01,
    -4.9486764754019213e-01,      7.0047448219202502e-01,      5.1423876868226370e-01,
    -6.2438983144173943e-01,      5.4170305683326148e-01,      5.6275672951076872e-01,
    -7.3428450139077006e-01,      3.8166488073486327e-01,      5.6138951702979889e-01,
    -8.0989114114122063e-01,      2.0921745876183720e-01,      5.4800035989970874e-01,
    -8.4928298013746839e-01,      1.6038410789326149e-02,      5.2769421924839510e-01,
    -7.9462472321999711e-01,     -1.7681271087024275e-01,      5.8078293236135903e-01,
    -6.9492222339230125e-01,     -3.7133740528338893e-01,      6.1578538052872012e-01,
    -5.8904908739445072e-01,     -5.5730626535313987e-01,      5.8517595579270032e-01,
    -4.7448547691242760e-01,     -6.9327993134680055e-01,      5.4242646413220086e-01,
    -3.3054165919756917e-01,     -8.0438332345888619e-01,      4.9366960659549908e-01,
    -1.6630776465645253e-01,     -8.7315269652004146e-01,      4.5819875160759044e-01,
     1.3405851065577399e-02,     -8.8666453160965697e-01,      4.6221887839273218e-01,
     1.6705677109295253e-01,     -8.4454662381276113e-01,      5.0875636157050985e-01,
     3.3720313237652355e-01,     -7.9391495649757848e-01,      5.0595759640992533e-01,
     5.4142548385441391e-01,     -7.0971727707195853e-01,      4.5073255047597649e-01,
     6.9685612784610396e-01,     -5.8073968622569350e-01,      4.2087165968477563e-01,
     7.8339540098176796e-01,     -4.2069816586839043e-01,      4.5749830486635390e-01,
     8.5674215792316377e-01,     -2.4008119953628232e-01,      4.5645798543390659e-01,
     8.9920050588696299e-01,     -3.5599185605149505e-02,      4.3608617060952509e-01,
     9.0560128560976216e-01,      1.4805360239690415e-01,      3.9744992429140719e-01,
     8.7375996501862463e-01,      3.1473194163680884e-01,      3.7079283763872278e-01,
     7.9579619301413962e-01,      4.6762304736000904e-01,      3.8475590282923655e-01,
     6.8891081640972596e-01,      6.1991902543347366e-01,      3.7563584618523410e-01,
     5.6952695117091179e-01,      7.4823724705526540e-01,      3.4026471167184469e-01,
     4.1544911613058222e-01,      8.3166364459313580e-01,      3.6842585979865800e-01,
     2.4172417072845084e-01,      9.0921286933187462e-01,      3.3897106591410148e-01,
     6.7624320110661149e-02,      9.5155397455868984e-01,      2.9995330441775642e-01,
    -1.1314049509153280e-01,      9.4205357627909447e-01,      3.1580735867014093e-01,
    -2.9838280569688963e-01,      9.0252499318404367e-01,      3.1050980329547451e-01,
    -4.9670722501559100e-01,      7.9849202058903956e-01,      3.4013589295007352e-01,
    -6.5914633352582441e-01,      6.3326075673623927e-01,      4.0559453272610585e-01,
    -7.8620962724541443e-01,      4.7603239119189827e-01,      3.9404007989384793e-01,
    -8.7051496734977962e-01,      3.2087212288008199e-01,      3.7315515858479248e-01,
    -9.1627448335879513e-01,      1.3436484462688189e-01,      3.7734223150075163e-01,
    -9.1582146868119763e-01,     -8.9918131073233265e-02,      3.9138953384947711e-01,
    -8.5049597846110603e-01,     -2.8185241873170658e-01,      4.4408986103779996e-01,
    -7.6846300034607651e-01,     -4.4761310238638413e-01,      4.5728232818592773e-01,
    -6.8355700077305415e-01,     -5.9855927407257981e-01,      4.1770398862813507e-01,
    -5.5548953676043067e-01,     -7.4288632633869656e-01,      3.7356563103245682e-01,
    -4.0289547912761525e-01,     -8.5651051959044933e-01,      3.2259101464459267e-01,
    -2.2900374310426258e-01,     -9.3358691432849761e-01,      2.7563156393786159e-01,
    -4.1994548417450282e-02,     -9.6188704453775953e-01,      2.7020320400326253e-01,
     1.4135197404517486e-01,     -9.4838408701592736e-01,      2.8387892300855661e-01,
     2.9480663752965230e-01,     -8.9989574746463374e-01,      3.2136690894292219e-01,
     4.6474211548098121e-01,     -8.2149464075698286e-01,      3.3039570412736569e-01,
     6.4402057280238423e-01,     -7.1788962459689487e-01,      2.6433310179282937e-01,
     7.8175301223319305e-01,     -5.7494296656414368e-01,      2.4145975454048463e-01,
     8.6762455772581071e-01,     -4.0800018517210118e-01,      2.8418915484342899e-01,
     9.3374257991234111e-01,     -2.0229988586328670e-01,      2.9529570033840036e-01,
     9.6945606925306804e-01,     -4.2872048741527122e-03,      2.4522754670460156e-01,
     9.6714860500721145e-01,      1.7196288641349938e-01,      1.8722270569816798e-01,
     9.2743067545239033e-01,      3.3744777651766872e-01,      1.6124931117125599e-01,
     8.5079457342757525e-01,      4.9298180174538825e-01,      1.8198224356804979e-01,
     7.5204832636565111e-01,      6.3397573499538051e-01,      1.8027224480682019e-01,
     6.4062114674999038e-01,      7.5646435801802303e-01,      1.3178095987322133e-01,
     4.9426815520496481e-01,      8.5479193810762333e-01,      1.5820788000758348e-01,
     3.3231874771470477e-01,      9.2817135635909342e-01,      1.6751770996480100e-01,
     1.5284763175113864e-01,      9.8159678610500622e-01,      1.1447860488489149e-01,
    -3.8568700996102045e-02,      9.9304378034726282e-01,      1.1124974434618060e-01,
    -2.2244733206529779e-01,      9.6670409199181162e-01,      1.2649262027216457e-01,
    -4.2655076662417329e-01,      8.8879872428732332e-01,      1.6760450828529830e-01,
    -6.2297807091470958e-01,      7.4695013048518633e-01,      2.3230115309130520e-01,
    -7.6247516534317461e-01,      6.0612072235966596e-01,      2.2638306509342804e-01,
    -8.6233172564433269e-01,      4.7404394210675188e-01,      1.7795037482162743e-01,
    -9.3387270359378971e-01,      3.0402771144417046e-01,      1.8827884681091006e-01,
    -9.7304472679398657e-01,      1.1024601141824235e-01,      2.0255808210186263e-01,
    -9.7316229579995184e-01,     -7.3066974862013134e-02,      2.1821173941353603e-01,
    -9.2706600410320039e-01,     -2.6991179631316153e-01,      2.6018502310303560e-01,
    -8.6211346477977802e-01,     -4.3689185496028254e-01,      2.5668245151308999e-01,
    -7.8331592779451675e-01,     -5.7687711943074516e-01,      2.3157924419235285e-01,
    -6.6040196234274184e-01,     -7.1826661838658012e-01,      2.1900299779994525e-01,
    -5.0549357418307705e-01,     -8.4687922265553817e-01,      1.6514184416424924e-01,
    -3.2083179560891845e-01,     -9.4008740121072965e-01,      1.1533705827362668e-01,
    -1.0121665458081842e-01,     -9.9101607981573536e-01,      8.7420354506940842e-02,
     9.9715545922655746e-02,     -9.9160736292610530e-01,      8.2290021825749002e-02,
     2.7720099421886329e-01,     -9.5656432510089917e-01,      9.0245779670492410e-02,
     4.3019139654600175e-01,     -8.9214454227561890e-01,      1.3788936879116984e-01,
     5.8336154432721576e-01,     -8.0404234400130115e-01,      1.1491395760767763e-01,
     7.3760983316028672e-01,     -6.7204307099766325e-01,      6.5496906410026035e-02,
     8.5771309983795774e-01,     -5.1031095921146907e-01,      6.2537694833051605e-02,
     9.3219480897613771e-01,     -3.3913609278961210e-01,      1.2648932241631206e-01,
     9.8411982875726733e-01,     -1.4312833248183399e-01,      1.0498782352128497e-01,
     9.9867839704839367e-01,      4.4927816435633455e-02,      2.4958977927330438e-02,
     9.7406488071920594e-01,      2.2309738794818162e-01,     -3.7751339581236637e-02,
     9.1761587516335152e-01,      3.9491997992534761e-01,     -4.4936790094083259e-02,
     8.3566975870930471e-01,      5.4853810953949744e-01,     -2.7604288825589785e-02,
     7.3270061569680200e-01,      6.7821136454961051e-01,     -5.6383976032935941e-02,
     5.9843270421532690e-01,      7.9766301349635249e-01,     -7.4913386156606102e-02,
     4.4921303769603727e-01,      8.9253510873683695e-01,     -3.9858831342893544e-02,
     2.7742724736313229e-01,      9.5983577778467555e-01,     -4.1825854504149509e-02,
     7.4355822578769201e-02,      9.9395512684206422e-01,     -8.0773866274993569e-02,
    -1.1568367585648862e-01,      9.8995970167701597e-01,     -8.1222387282598021e-02,
    -3.1801861899719081e-01,      9.4769137791314351e-01,     -2.7298538426204340e-02,
    -5.2315258953648225e-01,      8.5139185010915142e-01,      3.7990599218611754e-02,
    -6.7685669110187408e-01,      7.3472105657790132e-01,      4.5276801252677100e-02,
    -7.8646849828533949e-01,      6.1762905065456852e-01,     -1.2872421447274436e-03,
    -8.8316179570291220e-01,      4.6832938601462826e-01,     -2.6321641399594545e-02,
    -9.5608095103415314e-01,      2.9310267404535306e-01,      1.9374491713263148e-04,
    -9.9416742503490507e-01,      1.0750531593049320e-01,      8.5870860104812728e-03,
    -9.9744850636954363e-01,     -7.1343883263568802e-02,      2.5548898294941646e-03,
    -9.6986293019825343e-01,     -2.3463199175087951e-01,      6.5678954576584889e-02,
    -9.1183638659131505e-01,     -4.0765547758655135e-01,      4.8697183509674197e-02,
    -8.3811332974469088e-01,     -5.4549049803572014e-01,      2.4825505048482683e-03,
    -7.3523549584770376e-01,     -6.7701977675162284e-01,      3.2756488407093151e-02,
    -5.9396770091268714e-01,     -8.0442814575866228e-01,      9.8857768423679937e-03,
    -4.2189675920970193e-01,     -9.0579767907322650e-01,     -3.9162343570157597e-02,
    -2.1183067921902440e-01,     -9.7559486122118855e-01,     -5.7813753557578319e-02,
     1.3761537527447751e-02,     -9.9533260725399719e-01,     -9.5517647698427438e-02,
     2.0795493126984180e-01,     -9.7092540270736349e-01,     -1.1856900496377354e-01,
     3.8639700925520604e-01,     -9.1761988950291806e-01,     -9.3118685704237786e-02,
     5.3426685486947756e-01,     -8.4108468001655035e-01,     -8.4471822694516488e-02,
     6.7444988705480680e-01,     -7.3149332945642609e-01,     -1.0017414243461695e-01,
     8.0896552902772423e-01,     -5.7478098398168831e-01,     -1.2329474156644556e-01,
     9.1209482703407796e-01,     -3.9857545383638543e-01,     -9.6024132887492050e-02,
     9.7007558267049465e-01,     -2.3405230132612864e-01,     -6.4597864906229374e-02,
     9.9057918076926388e-01,     -4.8584740574644208e-02,     -1.2803284582398555e-01,
     9.6611719027890952e-01,      1.5097136643497555e-01,     -2.0934474238523776e-01,
     9.0537911034396812e-01,      3.4813899878043519e-01,     -2.4308003636851960e-01,
     8.2338107995986742e-01,      5.1322577610192488e-01,     -2.4216296147160740e-01,
     7.1226182106993308e-01,      6.5075463297373926e-01,     -2.6306179104795258e-01,
     5.7177762834491763e-01,      7.7472752037075721e-01,     -2.6994001723427191e-01,
     4.1183656872245494e-01,      8.8036854990599611e-01,     -2.3524850052514307e-01,
     2.2584908943740678e-01,      9.4651149145532465e-01,     -2.3045213243385509e-01,
     3.4823564914827160e-02,      9.6152178502195829e-01,     -2.7251270842808384e-01,
    -1.8422605218603977e-01,      9.5371056089847961e-01,     -2.3769082381668394e-01,
    -4.0009571251645604e-01,      9.0191671616134350e-01,     -1.6269498435627316e-01,
    -5.7295351396371841e-01,      8.0628791074964401e-01,     -1.4705126934372678e-01,
    -6.8944859386803281e-01,      6.9631084910555074e-01,     -1.9952904007010938e-01,
    -7.9433487653100210e-01,      5.6703776096254355e-01,     -2.1794559313980838e-01,
    -8.9443093722005607e-01,      3.9558847083819443e-01,     -2.0857387248538869e-01,
    -9.5832935520455520e-01,      2.1522698672244109e-01,     -1.8783554280167458e-01,
    -9.7941267816960065e-01,      3.2267305049844715e-02,     -1.9927274491377533e-01,
    -9.6972958106215956e-01,     -1.5448635517104281e-01,     -1.8909919534195554e-01,
    -9.3900142303604228e-01,     -3.1383816193572589e-01,     -1.4064826927158550e-01,
    -8.6538486734774611e-01,     -4.5378291881139782e-01,     -2.1257961793298080e-01,
    -7.7663995737010671e-01,     -5.9484029861626675e-01,     -2.0735331142344837e-01,
    -6.6866937995801978e-01,     -7.2365345565316719e-01,     -1.7090036988780014e-01,
    -5.2302383690740117e-01,     -8.2776730646834074e-01,     -2.0309444199388710e-01,
    -3.4373265928938246e-01,     -9.1310400100989342e-01,     -2.1929191110839863e-01,
    -1.4020413132158507e-01,     -9.6187344534200114e-01,     -2.3482392703101632e-01,
     5.4089549724683474e-02,     -9.5651335368418811e-01,     -2.8662959518236819e-01,
     2.5693921757208288e-01,     -9.1948784851251886e-01,     -2.9753039325632757e-01,
     4.2818864501135007e-01,     -8.5230965248851942e-01,     -3.0037100485606699e-01,
     5.6866049754066439e-01,     -7.6467622599661045e-01,     -3.0314271875205057e-01,
     7.0331442335207961e-01,     -6.4816259160844614e-01,     -2.9194875705225121e-01,
     8.2083071439201294e-01,     -4.8770516963773641e-01,     -2.9728875831306573e-01,
     9.0245074832657723e-01,     -3.2246465860157630e-01,     -2.8565572075099332e-01,
     9.4675154522186111e-01,     -1.5539248924146387e-01,     -2.8198348516776595e-01,
     9.3942350492749549e-01,      3.2543803956647688e-02,     -3.4121016868459697e-01,
     8.8668212713719485e-01,      2.4054756203022190e-01,     -3.9488185045248136e-01,
     7.9324461358276110e-01,      4.3784207415117576e-01,     -4.2315162899950209e-01,
     6.6487436468633632e-01,      6.0367471555904895e-01,     -4.3990785055241055e-01,
     5.0004863350302198e-01,      7.4806831645294913e-01,     -4.3628563814433624e-01,
     3.2757498934520318e-01,      8.5466883903562485e-01,     -4.0278505674489390e-01,
     1.4294731737798708e-01,      8.9093273557771968e-01,     -4.3105095421589995e-01,
    -7.6698562988817831e-02,      9.0652431185639948e-01,     -4.1512769414811601e-01,
    -2.9885304459404555e-01,      8.8587012659112152e-01,     -3.5485345785309719e-01,
    -4.7723375184225059e-01,      8.1165448767945980e-01,     -3.3684557104460039e-01,
    -5.9070112499490679e-01,      7.0029651018547101e-01,     -4.0082038215615007e-01,
    -7.0213088631686549e-01,      5.7574607854625981e-01,     -4.1896142008363613e-01,
    -8.1658626429664816e-01,      4.2840126654835037e-01,     -3.8685814943182228e-01,
    -8.9528114267496217e-01,      2.4000526116685136e-01,     -3.7532539240350543e-01,
    -9.2343505506915269e-01,      4.3346961980239497e-02,     -3.8129875420267828e-01,
    -9.1221413535555362e-01,     -1.5303667651190850e-01,     -3.8005939917295717e-01,
    -8.7018247282753725e-01,     -3.2052399576051893e-01,     -3.7422831550467073e-01,
    -7.7754566066867892e-01,     -4.5365801641942238e-01,     -4.3545051350722125e-01,
    -6.9307513323376513e-01,     -5.9990085243931679e-01,     -3.9970717648746290e-01,
    -5.7701457626285246e-01,     -7.2209212645182774e-01,     -3.8162434368954823e-01,
    -4.1846692182042255e-01,     -8.1743964713322326e-01,     -3.9582554065756287e-01,
    -2.3368486733364080e-01,     -8.8964388970066599e-01,     -3.9233293552482690e-01,
    -3.2625495764769953e-02,     -8.9680965749447417e-01,     -4.4121198448222915e-01,
     1.6947284787875871e-01,     -8.6927534635012482e-01,     -4.6436981605153183e-01,
     3.4593676955113895e-01,     -8.0092836208547513e-01,     -4.8871413963543159e-01,
     4.9103213011890401e-01,     -7.0842499523665725e-01,     -5.0697285264088232e-01,
     6.3193718428873569e-01,     -6.0825429379031171e-01,     -4.8029377384979055e-01,
     7.4281163437851261e-01,     -4.6966285248639178e-01,     -4.7712438716362920e-01,
     8.2190982368441334e-01,     -3.0601915800948737e-01,     -4.8043367561217115e-01,
     8.6616092257174337e-01,     -1.4064864862835849e-01,     -4.7956565123941358e-01,
     8.5553849536203119e-01,      3.7528742797396118e-02,     -5.1637726171639153e-01,
     8.0448878602670948e-01,      2.2480485137764658e-01,     -5.4978229505354703e-01,
     7.1145579332727904e-01,      4.1242549510756199e-01,     -5.6897791268759645e-01,
     5.6453522080296559e-01,      5.8830478759172433e-01,     -5.7896240065275661e-01,
     3.8670767623230412e-01,      7.1867321551973584e-01,     -5.7789789966527372e-01,
     2.1422842233063355e-01,      7.8406081520282322e-01,     -5.8254169046447823e-01,
     1.0361246385254807e-02,      8.1720823245021823e-01,     -5.7624938124820069e-01,
    -1.9615405632085922e-01,      8.1989326170776311e-01,     -5.3786487670703931e-01,
    -3.6938638067775181e-01,      7.7023869248417431e-01,     -5.1989042919644191e-01,
    -4.8760735565841418e-01,      6.5201177350003270e-01,     -5.8062011154036919e-01,
    -6.0237159364102477e-01,      5.2191612787279040e-01,     -6.0394703297610908e-01,
    -7.3010673652552094e-01,      3.9385167329445642e-01,     -5.5841294104203054e-01,
    -8.1880966713748538e-01,      2.1610187579389947e-01,     -5.3183710690450914e-01,
    -8.4031493486836450e-01,      1.0671763732338735e-03,     -5.4209747405034525e-01,
    -8.0854389266710569e-01,     -1.9454172138961801e-01,     -5.5534700167506823e-01,
    -7.2110482448365643e-01,     -3.4772295536683362e-01,     -5.9924667576662471e-01,
    -6.1839433346137351e-01,     -5.1389814750879614e-01,     -5.9455625665692158e-01,
    -5.0937990644362607e-01,     -6.4751421941261555e-01,     -5.6679577148206917e-01,
    -3.4951605400528346e-01,     -7.5195666242557169e-01,     -5.5892728134022096e-01,
    -1.6019485528309266e-01,     -8.1245139918663822e-01,     -5.6058927237372536e-01,
     6.3020179170373641e-02,     -7.9192147329951745e-01,     -6.0736219601194785e-01,
     2.4312793623186985e-01,     -7.2488018098470186e-01,     -6.4454443589190791e-01,
     3.9695485919687218e-01,     -6.2205801373103198e-01,     -6.7488566980859399e-01,
     5.4419171850857473e-01,     -5.2169282729226063e-01,     -6.5703269892638649e-01,
     6.5405822499545307e-01,     -3.8911401627379344e-01,     -6.4868954104030019e-01,
     7.2239826857434553e-01,     -2.2736848835706835e-01,     -6.5302703777332938e-01,
     7.4431715499691442e-01,     -5.7066157603096697e-02,     -6.6538366859558395e-01,
     7.1183483127693559e-01,      1.1902197699989876e-01,     -6.9218851620925670e-01,
     6.4353407858687206e-01,      2.8945375665394801e-01,     -7.0857632789718761e-01,
     5.4390246514374496e-01,      4.4629707466330587e-01,     -7.1062580135928988e-01,
     3.9747236621460202e-01,      5.7316938825363495e-01,     -7.1658395911764527e-01,
     2.2846945967430712e-01,      6.3676065459555720e-01,     -7.3643572343767216e-01,
     5.4820351142376905e-02,      6.9442605704546756e-01,     -7.1747277327917569e-01,
    -1.4181638180506936e-01,      7.0444271002330350e-01,     -6.9544847555138301e-01,
    -3.0601732280143301e-01,      6.4882792821099644e-01,     -6.9668911123891519e-01,
    -4.1368170488440204e-01,      5.1505486374488629e-01,     -7.5072360717954723e-01,
    -5.5496613439576581e-01,      3.8643120220367166e-01,     -7.3667056113112450e-01,
    -6.8404125222204615e-01,      2.4539267930793038e-01,     -6.8692794250967126e-01,
    -7.3620707946713104e-01,      5.9712019146120399e-02,     -6.7411691190176404e-01,
    -6.9788230645933258e-01,     -1.5200920950627697e-01,     -6.9989533971607432e-01,
    -5.7081295425104095e-01,     -3.7515175294802339e-01,     -7.3036547941371410e-01,
    -4.2175584030455604e-01,     -5.5753245093161474e-01,     -7.1503816494449057e-01,
    -2.5437101628435821e-01,     -6.6494288024768133e-01,     -7.0224379818007687e-01,
    -6.1745045543132725e-02,     -7.0557464965353789e-01,     -7.0594048128518860e-01,
     1.4127818461164954e-01,     -6.2661211196415534e-01,     -7.6641877305599548e-01,
     2.9905920466496561e-01,     -5.1285471896900281e-01,     -8.0470095646543338e-01,
     4.4897828834267994e-01,     -3.9928275611321556e-01,     -7.9936961242438509e-01,
     5.4758959028828991e-01,     -2.5842974487292436e-01,     -7.9583899601164187e-01,
     5.8777258638236207e-01,     -8.9065018383828243e-02,     -8.0410870483888996e-01,
     5.5838808643049387e-01,      8.3487607368053982e-02,     -8.2536813868021874e-01,
     4.8556321198342695e-01,      2.4631830462329174e-01,     -8.3878224824792569e-01,
     3.8137770834807788e-01,      3.9154156763609377e-01,     -8.3740446881315278e-01,
     2.2583761058221682e-01,      4.6484794361855913e-01,     -8.5610382720795719e-01,
     5.1491808076492389e-02,      5.2879109040862504e-01,     -8.4718863094677510e-01,
    -1.1872047261758756e-01,      5.5254658506375842e-01,     -8.2498346693484437e-01,
    -2.6215470209244257e-01,      4.5394929431735404e-01,     -8.5158966078716514e-01,
    -4.1111259088296742e-01,      3.1079648085824696e-01,     -8.5696673512081156e-01,
    -5.5295524558263875e-01,      1.7773596970243491e-01,     -8.1403342772676035e-01,
    -6.0098424742380618e-01,     -1.0538122123922494e-02,     -7.9919139280308982e-01,
    -5.4795660404932911e-01,     -2.1724110190821821e-01,     -8.0780558534862201e-01,
    -3.9939781639073091e-01,     -4.1228096124183322e-01,     -8.1884418130669157e-01,
    -2.1746581407024679e-01,     -5.2499370826492353e-01,     -8.2285492402549865e-01,
    -2.5970189173202155e-02,     -5.5722924522351891e-01,     -8.2995247908656522e-01,
     1.2560622645929923e-01,     -4.4149504272983514e-01,     -8.8842850197392720e-01,
     2.8491645095032458e-01,     -3.2047950105934336e-01,     -9.0339111429027552e-01,
     3.8911999460025343e-01,     -1.6753114469052957e-01,     -9.0582500813401012e-01,
     3.9333108248868054e-01,      1.8469652535215329e-03,     -9.1939504472649569e-01,
     3.1957966926054959e-01,      1.7791345374272430e-01,     -9.3070706345909571e-01,
     1.8293770781829197e-01,      2.8092384011411042e-01,     -9.4213353146660206e-01,
     1.4497531695760716e-02,      3.5058809443447775e-01,     -9.3641754021138035e-01,
    -1.4901751757979714e-01,      3.2382842934670192e-01,     -9.3430665619014130e-01,
    -3.0333495673084798e-01,      1.9570498700964184e-01,     -9.3257035235131214e-01,
    -4.3335585511223695e-01,      4.4061699154812436e-02,     -9.0014513802360374e-01,
    -4.0927666999556650e-01,     -1.4440065510421615e-01,     -9.0091123769371051e-01,
    -3.0218101545001530e-01,     -2.8956186858946509e-01,     -9.0820733214424942e-01,
    -1.3157193529018418e-01,     -3.8450470804206721e-01,     -9.1369850352152826e-01,
     2.7418992969621237e-02,     -3.0413798346172283e-01,     -9.5223331481332274e-01,
     1.7664475615462172e-01,     -1.7261969133641308e-01,     -9.6901964494327752e-01,
     1.9768573722830793e-01,      1.2908281633725205e-02,     -9.8018045561098843e-01,
     5.0916808773057765e-02,      1.2770583017479939e-01,     -9.9050426527286251e-01,
    -1.1623817515301313e-01,      1.1813763474578851e-01,     -9.8617046492671245e-01,
    -2.4998426354769901e-01,     -1.4523416772369846e-02,     -9.6814097028468460e-01,
    -1.5293200167983850e-01,     -1.6844979092376627e-01,     -9.7377434285358810e-01,
    -1.4275241161755283e-03,     -8.4668893443673982e-02,     -9.9640811952629205e-01
};


static constexpr auto weights = detail::create_array<339, T>(4.0 * M_PI / 339.0);
};
}  // namespace WomersleyGrids
}  // namespace IntegratorXX
