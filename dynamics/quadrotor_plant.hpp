#pragma once

#include <vector>
#include <casadi/casadi.hpp>
#include "roll_pitch_yaw.hpp"

namespace quadrotor {

using casadi::SX;

class DynamicSystem {

  public:
    SX&                     X()       { return X_; }
    [[nodiscard]] const SX& X() const { return X_; }

    SX&                     U()       { return U_; }
    [[nodiscard]] const SX& U() const { return U_; }

    SX&                     ode()       { return ode_; }
    [[nodiscard]] const SX& ode() const { return ode_; }

    [[nodiscard]] int nx() const { return static_cast<int>(X_.numel()); }
    [[nodiscard]] int nu() const { return static_cast<int>(U_.numel()); }

  protected:
    virtual void set_ode() {};

  private:
    SX X_;   // state vector
    SX U_;   // control vector
    SX ode_; // dX/dt = f(x, u)
};

class Quadrotor : public DynamicSystem {

  public:
    Quadrotor();
    Quadrotor(double m, double L, double kF, double kM, SX I_);

    double&                     g()       { return g_; }
    [[nodiscard]] const double& g() const { return g_; }

    double&                     m()  { return m_; }
    [[nodiscard]] const double& m() const { return m_; }

    double&                     L()       { return L_; }
    [[nodiscard]] const double& L() const { return L_; }

    double&                     kF()       { return kF_; }
    [[nodiscard]] const double& kF() const { return kF_; }

    double&                     kM()       { return kM_; }
    [[nodiscard]] const double& kM() const { return kM_; }

    SX&                     I()       { return I_; }
    [[nodiscard]] const SX& I() const { return I_; }

  protected:
    void set_ode() override;

  private:
    double g_;  // gravitational acceleration (m/s^2)
    double m_;  // mass (kg)
    double L_;  // length of arms (m)
    double kF_; // force input constant
    double kM_; // moment input constant
    SX I_;      // moment of intertia about center of mass
};

//class PlanarQuadrotor : public Quadrotor {
//
//  public:
//    PlanarQuadrotor();
//    PlanarQuadrotor(double m_arg, double l_arg);
//
//    SX calculate_ode() override;
//};
//
//class Crazyflie : Quadrotor {
//  // crazyflie 2.1 specs
//  public:
//    Crazyflie();
//    Crazyflie(double m_arg, double l_arg);
//
//    SX calculate_ode() override;
//};
//
//class PlanarCrazyflie : Crazyflie {
//  // crazyflie 2.1 specs
//  public:
//    PlanarCrazyflie();
//    PlanarCrazyflie(double m_arg, double l_arg);
//
//    SX calculate_ode() override;
//};

} // namespace quadrotor

