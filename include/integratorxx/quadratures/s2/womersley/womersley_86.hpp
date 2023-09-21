#pragma once

namespace IntegratorXX {
namespace WomersleyGrids {

/**
 *  \brief Womersley Quadrature specification for index 12 grid with 86 points
 * 
 */
template <typename T>
struct womersley_86 {

  static constexpr std::array<cartesian_pt_t<T>,86> points = {
     0.0000000000000000e+00,      0.0000000000000000e+00,      1.0000000000000000e+00,
     9.5516280663209729e-01,      0.0000000000000000e+00,     -2.9608109163993362e-01,
     9.4207955549233324e-01,      3.3009429934076029e-01,     -5.9362148428946986e-02,
    -8.0438154482809987e-01,      5.2810115812758374e-01,     -2.7217548957293042e-01,
     3.3567061410450411e-01,      5.9794905084362837e-01,      7.2786136827139625e-01,
    -8.8832210433867276e-01,      4.3453421838484185e-01,      1.4853905882287913e-01,
     7.4260463414904165e-01,      1.8112230673114985e-01,     -6.4477365590163127e-01,
     3.2917595739443783e-01,     -2.2407163711111855e-01,      9.1729771095092050e-01,
    -1.8650239433050182e-02,     -7.7904672910476547e-01,     -6.2668840937124104e-01,
     6.2519804309412330e-01,      6.7338955867630055e-03,      7.8043709647959802e-01,
    -9.6121753193579032e-02,      7.7394490885945555e-01,      6.2591523916064051e-01,
    -4.8670179721916262e-01,     -8.1710973271716192e-01,     -3.0895476251795390e-01,
    -7.9117860273360885e-01,     -5.0856874958206198e-02,     -6.0946697765012658e-01,
     1.5819165083650261e-01,      3.9299907545783253e-01,     -9.0582952496311964e-01,
    -6.4410742931972553e-01,      1.8296936136465711e-01,      7.4272998613018459e-01,
    -6.6625730070885669e-01,      5.1271068143836429e-01,      5.4150620161837137e-01,
     8.8430662313444242e-01,     -4.4544721060059056e-01,     -1.3992347497368893e-01,
    -3.8908160834754840e-01,     -6.7068798420346487e-01,     -6.3150069666689823e-01,
     3.3697039742509044e-01,     -9.1502854501080844e-01,      2.2175146690513606e-01,
     5.1469289348355873e-01,     -1.7639543410447664e-01,     -8.3903270271462915e-01,
    -8.9607045535351049e-01,     -4.4379173202483641e-01,      1.0326549711687219e-02,
     7.0678408796600134e-01,     -3.3478700184495525e-01,      6.2319653111497142e-01,
    -8.9103881062487689e-01,     -2.8001390999047621e-01,      3.5727027328347682e-01,
     6.6905253966941969e-01,      6.0923444586895348e-01,      4.2567838696444155e-01,
    -4.7030978850856725e-01,      2.2814209113628964e-01,     -8.5250213435802469e-01,
    -4.5272456436071057e-02,     -3.9150332357338818e-01,      9.1906232232598617e-01,
     1.2469156567170592e-01,     -1.4910312123512653e-01,     -9.8092827091907275e-01,
     2.1547589910863515e-01,     -8.0023721028169703e-01,      5.5963429503908391e-01,
     4.4263179691676591e-01,      8.8130914611573397e-01,      1.6544268291765174e-01,
    -1.0883233201205400e-01,     -9.4721722818045107e-01,     -3.0155438339868551e-01,
    -2.4420991262447247e-01,     -4.1611498201843305e-01,     -8.7590515486309739e-01,
     8.1267068216668881e-01,      4.4252395645745018e-01,     -3.7912914726775265e-01,
    -4.4695700783132192e-01,     -1.3443779486982377e-01,     -8.8439578948625186e-01,
    -3.3948271153722021e-01,     -9.3896962714342558e-01,      5.5565527708038416e-02,
     3.2410908664497917e-01,     -8.0145605121022434e-01,     -5.0261466147802691e-01,
    -1.8115177470688795e-01,      4.9396442203800356e-01,     -8.5040177815030826e-01,
    -9.8998144229801230e-02,      9.5042493107766035e-01,     -2.9477418106930564e-01,
    -6.5110236222632745e-01,     -7.5480652486473376e-01,      7.9579041995493219e-02,
    -5.0264944382434285e-01,     -4.8228066428614702e-01,      7.1746003197306507e-01,
    -1.3878908906763374e-01,      9.4292977295044667e-01,      3.0268966292127564e-01,
    -6.6454628468833032e-01,     -3.8840138015573622e-01,     -6.3837497084398276e-01,
    -6.7777254369133988e-01,     -6.1219079510467334e-01,      4.0724293659593247e-01,
    -7.4416958731804794e-01,      3.0323464410553608e-01,     -5.9519776202962593e-01,
    -9.4290689023707530e-01,     -1.7442524283366692e-01,     -2.8372950323479657e-01,
    -3.6136600494889798e-01,      9.3237938655603625e-01,     -9.1263351167203549e-03,
     8.8331859840263804e-01,     -3.8582577051765682e-01,      2.6624561690374443e-01,
     4.8137556795853248e-01,     -5.1183720111977282e-01,     -7.1154777922672041e-01,
     1.2389800107889963e-01,     -4.9975022448054734e-01,     -8.5726250265615578e-01,
    -4.9169774767836077e-01,      7.7369505139504025e-01,      3.9953634675064698e-01,
     6.0167107935854092e-01,      7.5254081374196191e-01,     -2.6772044359016417e-01,
    -3.1830926920317631e-01,      1.7831843882528875e-01,      9.3106484388266508e-01,
    -4.8720528558805054e-01,      6.1134098621042854e-01,     -6.2361302766565641e-01,
     9.9416827041354860e-01,     -9.7871228717060735e-02,      4.5284353725607657e-02,
    -9.1500891448587793e-02,     -6.7034333965383863e-01,      7.3638807285686814e-01,
     4.2388270558574859e-01,      1.7057165318317524e-01,     -8.8951040636726908e-01,
     7.2585330397840808e-01,     -5.4576268904428960e-01,     -4.1866462514854874e-01,
    -6.8879215092801882e-01,      7.2114641958289516e-01,      7.4251022502887926e-02,
     8.0892481063413568e-01,     -1.8964758212658520e-01,     -5.5648400276563903e-01,
     6.5226385872147385e-01,      3.5383732111226179e-01,      6.7033648923049927e-01,
     1.3579871288592421e-01,     -9.8736576240231178e-01,     -8.1655133422387186e-02,
    -1.0261388431521144e-01,      7.9067036317634687e-01,     -6.0358161630414831e-01,
    -7.6626146336963896e-01,     -5.5719552014144680e-01,     -3.1996331366099140e-01,
    -5.1237386864842149e-01,      7.9262553082947185e-01,     -3.3048114409682738e-01,
     1.1886224779557640e-01,      9.9212174514951701e-01,      3.9575356605592832e-02,
     5.3533537104770357e-01,      5.2994879567413034e-01,     -6.5770077882632649e-01,
     2.3557931304068261e-02,      4.4864432372725477e-01,      8.9339985149986845e-01,
     7.4802040688164673e-01,      6.6037366065861358e-01,      6.6121851130757345e-02,
     3.1545503972944389e-01,      2.1326491887406265e-01,      9.2466544884457236e-01,
    -9.9529609036279554e-01,     -1.3225954614801000e-02,      9.5972739009892782e-02,
     3.0966556313699239e-01,      9.0265899038191921e-01,     -2.9885445636587654e-01,
     8.8960578357380093e-01,     -8.6713880983188656e-03,      4.5664686231320145e-01,
     3.7868361167186193e-01,     -5.3241536018391267e-01,      7.5705522023917715e-01,
    -1.0361523041444862e-01,      7.8707530872592316e-02,     -9.9149836540969682e-01,
    -3.3684760212483705e-01,      5.3996388232677428e-01,      7.7134473403617976e-01,
    -7.4621851119147287e-01,     -1.5233268645255751e-01,      6.4803756541833946e-01,
    -8.9356616821260304e-01,      1.9150717197816608e-01,      4.0603510452517638e-01,
    -7.3257524912105382e-02,     -9.5192790349887935e-01,      2.9743335654223546e-01,
     9.0892522859731273e-01,      2.8693295920312278e-01,      3.0253000800294488e-01,
     2.6730980347132566e-01,      8.3615554833679451e-01,      4.7894610130338877e-01,
     2.5281309121179418e-01,      7.1410498354159890e-01,     -6.5279369895319128e-01,
     4.9469331204964989e-01,     -8.4302889206641307e-01,     -2.1114169212787823e-01,
    -9.5404066739523230e-01,      2.0335560082470686e-01,     -2.2012020481837269e-01,
     6.8599369224278017e-01,     -7.2426661023353112e-01,      6.9645757257337307e-02,
     6.1221896606149240e-01,     -6.6510462446779894e-01,      4.2757897060793959e-01,
    -3.9472137008217872e-01,     -2.0335065684463793e-01,      8.9601537395365161e-01,
    -3.4501808581940013e-01,     -7.9948430995452158e-01,      4.9171877998919211e-01
};


static constexpr auto weights = detail::create_array<86, T>(4.0 * M_PI / 86.0);
};
}  // namespace WomersleyGrids
}  // namespace IntegratorXX