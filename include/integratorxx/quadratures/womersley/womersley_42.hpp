#pragma once

namespace IntegratorXX {
namespace WomersleyGrids {

/**
 *  \brief Womersley Quadrature specification for index 8 grid with 42 points
 * 
 */
template <typename T>
struct womersley_42 {

  static constexpr std::array<cartesian_pt_t<T>,42> points = {
     0.0000000000000000e+00,      0.0000000000000000e+00,      1.0000000000000000e+00,
     9.3248216132655637e-01,      0.0000000000000000e+00,      3.6121602789432539e-01,
    -9.7956656131034692e-01,     -8.2814270168621987e-02,      1.8327888208700083e-01,
     3.9738932087753187e-01,      9.0173013794260004e-01,     -1.7018955895857293e-01,
     6.0220344466110000e-02,     -9.9785807088911793e-01,     -2.5549568955875734e-02,
    -5.6709982693717609e-01,     -8.2364885309003422e-01,     -5.9421485694216667e-04,
     4.2414717665656149e-01,      7.8819299514543406e-01,      4.4592709598984598e-01,
     3.3514055813887228e-01,     -7.9436923438049611e-01,      5.0661457318173309e-01,
    -6.3486255483510967e-01,      7.2420784815945216e-01,      2.6920722340326109e-01,
    -5.1144426259146991e-01,      6.8484805373262084e-02,     -8.5658309444866698e-01,
     7.8070285259453875e-01,      5.5103704624743666e-01,     -2.9472229235951969e-01,
    -6.0607443503881864e-01,     -4.4232708952261013e-01,     -6.6107527942499389e-01,
     4.5436938471651184e-01,      4.9031139597450757e-01,     -7.4373597278191950e-01,
    -7.2684950185476893e-01,     -5.2405360534112433e-01,      4.4391172588978023e-01,
     8.0538367492989960e-01,      5.3859751608412820e-01,      2.4752747690795279e-01,
     6.0529235093974110e-01,     -4.8335320355514827e-01,     -6.3244829868283836e-01,
    -8.8214066699744209e-01,     -4.2377167279997480e-01,     -2.0553688953961324e-01,
    -5.4935182340687749e-01,     -2.5170247526455547e-01,      7.9678004371673028e-01,
     5.2173215746075685e-01,     -8.4895680483704083e-01,     -8.4070799878593366e-02,
    -8.9079994946924845e-01,      1.0157707288321312e-01,     -4.4289676933802835e-01,
     7.5178914827606835e-01,     -5.0532758529938782e-01,      4.2362378128456574e-01,
    -1.1260009749008823e-01,     -3.8179903075102317e-01,     -9.1736073502347038e-01,
    -7.6115438426043805e-02,      3.7504026947793434e-01,     -9.2387836661712641e-01,
    -8.7765172836372618e-01,      4.6488290475910576e-01,     -1.1666759859918491e-01,
     1.9337720432514655e-01,      5.0009572934043223e-01,      8.4410278896758284e-01,
     6.2934353923671515e-01,      2.0622043898635703e-01,      7.4926620113633802e-01,
    -3.9295603363840315e-01,      8.9254549935326688e-01,     -2.2124214610105705e-01,
     5.3243714845248889e-02,      8.2364020993449738e-01,     -5.6460775004292574e-01,
    -2.3194023438401173e-01,      7.6061029634195243e-01,      6.0636268418537842e-01,
    -8.2522349799878125e-01,      2.4602482853694707e-01,      5.0840727974137168e-01,
     4.4913957826789591e-01,     -2.8151685821646655e-01,      8.4795158928636150e-01,
    -1.7129970050895518e-02,      9.8775935573368312e-01,      1.5504199201066343e-01,
    -4.3001618057787228e-01,      3.0424358181475220e-01,      8.5001289835257732e-01,
    -2.3025372842153186e-01,     -8.6885338124455680e-01,      4.3826592663345904e-01,
     7.3429958698090647e-01,      2.4109574854103540e-02,     -6.7839726190486971e-01,
     2.0168965492904883e-01,     -7.5659508418765864e-01,     -6.2200093382379340e-01,
     8.7531209671403842e-01,     -4.6473801038936263e-01,     -1.3359384359097765e-01,
    -2.9032439082635958e-01,     -8.4143896543788510e-01,     -4.5573261407773408e-01,
    -5.2553721839997014e-01,      5.8891161460766173e-01,     -6.1399816144400543e-01,
     9.8183172842281263e-01,      4.8599913762707996e-02,     -1.8342438617732823e-01,
     2.4999831162012423e-01,     -7.1330817010645461e-02,     -9.6561522291836366e-01,
    -7.8946444706143440e-02,     -5.4236802292772512e-01,      8.3642356887747393e-01
};


static constexpr auto weights = detail::create_array<42, T>(4.0 * M_PI / 42.0);
};
}  // namespace WomersleyGrids
}  // namespace IntegratorXX