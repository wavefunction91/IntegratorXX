#include "SlaterAtomicShell.hpp"
static const std::vector<SlaterTypeAtom> k99l_neutral{
    SlaterTypeAtom(
        1,
        std::vector<SlaterTypeAtomicShell>{
            /**
             *
             * Slater-type orbital Hartree-Fock calculation for H
             *
             * Etot=-0.5 Ekin=0.5
             */
            SlaterTypeAtomicShell(0, std::vector<double>{1.000000},
                                  std::vector<int>{1},
                                  std::vector<std::vector<double>>{{1.0000000}},
                                  std::vector<int>{1}, std::vector<int>{0})}),
};
