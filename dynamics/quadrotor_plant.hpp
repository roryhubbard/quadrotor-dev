#pragma once

#include <vector>
#include <casadi/casadi.hpp>
#include "roll_pitch_yaw.hpp"

namespace quadrotor {

using casadi::SX;

template <typename T>
struct ThreeDOF {
  ThreeDOF()
      : x(T::sym("x")),
        y(T::sym("y")),
        z(T::sym("z")),
        xdt(T::sym("xdt")),
        ydt(T::sym("ydt")),
        zdt(T::sym("zdt")) {}
  T x;
  T y;
  T z;
  T xdt;
  T ydt;
  T zdt;
  [[nodiscard]] T  positions() const { return vertcat(x, y, z); };
  [[nodiscard]] T velocities() const { return vertcat(xdt, ydt, zdt); };
  [[nodiscard]] T full_state() const { return vertcat(positions(),
                                                      velocities()); };
};

template <typename T>
struct SixDOF {
  SixDOF()
      : xyz(),
        rpy() {}
  [[nodiscard]] T  positions() const { return vertcat(xyz.positions(),
                                                      rpy.positions()); };
  [[nodiscard]] T velocities() const { return vertcat(xyz.velocities(),
                                                      rpy.velocities()); };
  [[nodiscard]] T full_state() const { return vertcat(xyz.positions(),
                                                      rpy.positions(),
                                                      xyz.velocities(),
                                                      rpy.velocities()); };
  ThreeDOF<T> xyz;
  ThreeDOF<T> rpy;
};

class QuadRotor {

  public:
    QuadRotor();
    QuadRotor(double m, double L, SX I, double kF, double kM);
    [[nodiscard]] SixDOF<SX> X() const { return X_; } 
    [[nodiscard]] SX U() const { return U_; }
    [[nodiscard]] SX ode() const { return ode_; }

  private:
    const double g_;  // gravitational acceleration (m/s^2)
    const double m_;  // mass (kg)
    const double L_;  // length of arms (m)
    const double kF_; // force input constant
    const double kM_; // moment input constant
    const SX I_;       // moment of intertia about center of mass

    SixDOF<SX> X_; // state vector
    SX U_; // control vector (rotational velocity squared of eadh propeller)
    const SX ode_; // dX/dt = f(x, u)

    SX calculate_ode();
};

class PlanarQuadrotor: public QuadRotor {

  public:
    PlanarQuadrotor();
    PlanarQuadrotor(double m_arg, double l_arg);
};

class CrazyFlie : QuadRotor {
  // crazyflie 2.1 specs
};

class PlanarCrazyFlie final : CrazyFlie {
  // crazyflie 2.1 specs
};

} // namespace quadrotor

