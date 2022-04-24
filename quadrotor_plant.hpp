#pragma once

#include <casadi/casadi.hpp>

namespace quadrotor {

inline constexpr double kGimbalLockToleranceCosPitchAngle = 0.008;

template <typename T>
struct ThreeDOF {
  T x;
  T y;
  T z;
  T xdt;
  T ydt;
  T zdt;
  T positions() const { return T::sym({x, y, z}); };
  T velocities() const { return T::sym({xdt, ydt, zdt}); };
  T vertcat() const { return T::sym({positions(),
                                     velocities()}); };
};

template <typename T>
struct SixDOF {
  T vertcat() const { return T::sym({xyz.positions(),
                                     rpy.positions(),
                                     xyz.velocities(),
                                     rpy.velocities()}); };
  ThreeDOF<T> xyz;
  ThreeDOF<T> rpy;
};

template <typename T>
class QuadRotor {

  public:
    QuadRotor();
    QuadRotor(double m, double L, T I, double kF, double kM);
    T X() const { return X_; }
    T U() const { return U_; }
    T ode() const { return ode_; }

  private:
    const double g_;  // gravitational acceleration (m/s^2)
    const double m_;  // mass (kg)
    const double L_;  // length of arms (m)
    const double kF_; // force input constant
    const double kM_; // moment input constant
    const T I_;       // moment of intertia about center of mass

    SixDOF<T> X_; // state vector
    T U_; // control vector (rotational velocity squared of eadh propeller)
    T ode_; // dX/dt = f(x, u)

    T calculate_ode();
};

template <typename T>
class PlanarQuadrotor: public QuadRotor<T> {

  public:
    PlanarQuadrotor();
    PlanarQuadrotor(double m_arg, double l_arg);
};

template <typename T>
class CrazyFlie final : QuadRotor<T> {
  // crazyflie 2.1 specs
};

template <typename T>
class PlanarCrazyFlie final : CrazyFlie<T> {
  // crazyflie 2.1 specs
};

} // namespace quadrotor

