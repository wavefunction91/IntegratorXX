#include "quad_matcher.hpp"
#include <sstream>
#include <iomanip>

namespace IntegratorXX::Matchers {
    WithinAbsMatcher::WithinAbsMatcher(std::string msg, double target, double margin)
        :m_msg{msg}, m_target{ target }, m_margin{ margin } {
    }
    bool WithinAbsMatcher::match(double const& matchee) const {
        return (matchee + m_margin >= m_target) && (m_target + m_margin >= matchee);
    }

    std::string WithinAbsMatcher::describe() const {
        std::stringstream ss;
        ss << std::scientific;
        ss << std::setprecision(6);
        ss << " is not within " << m_margin << " of " << m_target << "\n(" << m_msg << ")";
        return ss.str();
    }

WithinAbsMatcher WithinAbs(std::string msg, double target, double margin) {
    return WithinAbsMatcher(msg, target, margin);
}
}
