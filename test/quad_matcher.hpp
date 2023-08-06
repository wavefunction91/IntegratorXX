#pragma once
#include <catch2/matchers/catch_matchers.hpp>

namespace IntegratorXX::Matchers {

class WithinAbsMatcher final : public Catch::Matchers::MatcherBase<double> {
  public:
    WithinAbsMatcher(std::string msg, double target, double margin);
    bool match(double const& matchee) const override;
    std::string describe() const override;
  public:
    std::string m_msg;
    double m_target;
    double m_margin;
};

WithinAbsMatcher WithinAbs(std::string, double, double);

}

