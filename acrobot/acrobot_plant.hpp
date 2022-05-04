#pragma once

#include <casadi/casadi.hpp>
#include "system.hpp"

namespace acrobot {

using casadi::SX;

class Acrobot : public systems::System {

  public:
    Acrobot();

    double&                     g()       { return g_; }
    [[nodiscard]] const double& g() const { return g_; }

    double&                     m1()  { return m1_; }
    [[nodiscard]] const double& m1() const { return m1_; }

    double&                     m2()  { return m2_; }
    [[nodiscard]] const double& m2() const { return m2_; }

    double&                     l1()       { return l1_; }
    [[nodiscard]] const double& l1() const { return l1_; }

    double&                     l2()       { return l2_; }
    [[nodiscard]] const double& l2() const { return l2_; }

    SX&                     I1()       { return I1_; }
    [[nodiscard]] const SX& I1() const { return I1_; }

    SX&                     I2()       { return I2_; }
    [[nodiscard]] const SX& I2() const { return I2_; }

  protected:
    void set_ode() override;

  private:
    double g_;  // gravitational acceleration (m/s^2)
    double m1_;  // link 1 mass (kg)
    double m2_;  // link 2 mass (kg)
    double l1_;  // link 1 length (m)
    double l2_;  // link 2 length (m)
    SX I1_;      // link 1 moment of intertia
    SX I2_;      // link 2 moment of intertia

};

} // namespace acrobot
