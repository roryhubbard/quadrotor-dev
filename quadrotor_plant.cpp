#include "quadrotor_plant.hpp"

using casadi::SX;

namespace quadrotor {

template <typename T>
void default_moment_of_inertia(T& I) {
  using casadi::Slice;
  Slice all;
  I(0, all) = {0.0015, 0, 0};
  I(1, all) = {0, 0.0025, 0};
  I(2, all) = {0, 0, 0.0035};
}

template <typename T>
T rotation_matrix(const T& rpy) {
  auto& sin = T::sin;
  auto& cos = T::sin;
  const auto& r = rpy(0);
  const auto& p = rpy(1);
  const auto& y = rpy(2);
  const auto cr = cos(r), cp = cos(p), cy = cos(y);
  const auto sr = sin(r), sp = sin(p), sy = sin(y);
  const auto Rxx = cy * cp;
  const auto Rxy = cy * sp * sr - sy * cr;
  const auto Rxz = cy * sp * cr + sy * sr;
  const auto Ryx = sy * cp;
  const auto Ryy = sy * sp * sr + cy * cr;
  const auto Ryz = sy * sp * cr - cy * sr;
  const auto Rzx = -sp;
  const auto Rzy = cp * sr;
  const auto Rzz = cp * cr;
  return T::sym({{Rxx, Rxy, Rxz},
                 {Ryx, Ryy, Ryz},
                 {Rzx, Rzy, Rzz}});
}

template <typename T>
T matrix_relating_angular_velocity_in_child_to_rpyDt(const T& rpy) {
  auto& sin = T::sin;
  auto& cos = T::cos;
  const auto& r = rpy(0);
  const auto& p = rpy(1);
  const auto sr = sin(r), cr = cos(r);
  const auto sp = sin(p), cp = cos(p);
  return T::sym({{1,   0,      -sp},
                 {0,  cr,  sr * cp},
                 {0, -sr,  cr * cp}});
}

template <typename T>
T angular_velocity_in_child_from_rpyDt(const T& rpy, const T& rpyDt) {
  const auto M = matrix_relating_angular_velocity_in_child_to_rpyDt(rpy);
  return M * rpyDt;
}

template <typename T>
T matrix_relating_rpyDt_to_angular_velocity_in_parent(const T& rpy) {
  auto& sin = T::sin;
  auto& cos = T::cos;
  const auto& p = rpy(1);
  const auto& y = rpy(2);
  const auto& sp = sin(p), cp = cos(p);
  if (std::abs(cp) < kGimbalLockToleranceCosPitchAngle) {
    throw std::domain_error("gimbal lock error");
  }
  const auto one_over_cp = 1 / cp;
  const auto sy = sin(y), cy = cos(y);
  const auto cy_over_cp = cy * one_over_cp;
  const auto sy_over_cp = sy * one_over_cp;
  return T::sym({{     cy_over_cp,      sy_over_cp, 0,
                              -sy,              cy, 0,
                  cy_over_cp * sp, sy_over_cp * sp, 1}});
}

template <typename T>
T matrix_relating_angular_velocity_in_parent_to_rpyDt(
    const T& rpy, const T& rpyDt) {
  auto& sin = T::sin;
  auto& cos = T::cos;
  const auto& p = rpy(1);
  const auto& y = rpy(2);
  const auto sp = sin(p), cp = cos(p);
  const auto sy = sin(y), cy = cos(y);
  const auto& pDt = rpyDt(1);
  const auto& yDt = rpyDt(2);
  const auto sp_pDt = sp * pDt;
  const auto cp_yDt = cp * yDt;
  return T::sym({{-cy * sp_pDt - sy * cp_yDt, -cy * yDt,  T(0),
                  -sy * sp_pDt + cy * cp_yDt, -sy * yDt,  T(0),
                                   -cp * pDt,       T(0), T(0)}});
}

template <typename T>
T rpyDDT_from_rpyDt_and_angular_accel_in_parent(
    const T& rpy, const T& rpyDt, const T& alpha_AD_A) {
  const T Minv = matrix_relating_rpyDt_to_angular_velocity_in_parent(rpy);
  const T MDt = matrix_relating_angular_velocity_in_parent_to_rpyDt(rpy, rpyDt);
  return Minv * (alpha_AD_A - MDt * rpyDt);
}

template <typename T>
QuadRotor<T>::QuadRotor()
    : QuadRotor(
        0.1, // m (kg)
        0.1, // L (m)
        default_moment_of_inertia(T::sym("I", 3, 3)),
        1.0, // kF
        0.1  // kM
        )
      {}

template <typename T>
QuadRotor<T>::QuadRotor(double m, double L, T I, double kF, double kM)
    : g_(9.81),
      m_(m),
      L_(L),
      I_(I),
      kF_(kF),
      kM_(kM),
      X_(SixDOF<T>()),
      U_(T::sym(4)),
      ode_(calculate_ode())
      { }

template <typename T>
T QuadRotor<T>::calculate_ode() {
  // https://andrew.gibiansky.com/downloads/pdf/Quadcopter%20Dynamics,%20Simulation,%20and%20Control.pdf
  // https://github.com/RobotLocomotion/drake/blob/master/examples/quadrotor/quadrotor_plant.cc#L62

  // For each rotor, calculate the Bz measure of its aerodynamic force on B.
  // Note: B is the quadrotor body and Bz is parallel to each rotor's spin axis.
  const auto uF_Bz = kF_ * U_;
  const auto thrust = T::sum1(uF_Bz);

  // Compute the net aerodynamic force on B (from the 4 rotors), expressed in B.
  const auto Faero_B = T::sym({0, 0, thrust});

  // Compute the Bx and By measures of the moment on B about Bcm (B's center of
  // mass) from the 4 rotor forces.  These moments arise from the cross product
  // of a position vector with an aerodynamic force at the center of each rotor.
  // For example, the moment of the aerodynamic forces on rotor 0 about Bcm
  // results from Cross( L_* Bx, uF_Bz(0) * Bz ) = -L_ * uF_Bz(0) * By.
  const auto Mx = L_ * (uF_Bz(1) - uF_Bz(3));
  const auto My = L_ * (uF_Bz(2) - uF_Bz(0));

  // For rotors 0 and 2, get the Bz measure of its aerodynamic torque on B.
  // For rotors 1 and 3, get the -Bz measure of its aerodynamic torque on B.
  // Sum the net Bz measure of the aerodynamic torque on B.
  // Note: Rotors 0 and 2 rotate one way and rotors 1 and 3 rotate the other.
  const auto uM_Bz = kM_ * U_;
  const auto Mz = uM_Bz(0) - uM_Bz(1) + uM_Bz(2) - uM_Bz(3);

  // For rotors 0 and 2, get the Bz measure of its aerodynamic torque on B.
  // For rotors 1 and 3, get the -Bz measure of its aerodynamic torque on B.
  // Sum the net Bz measure of the aerodynamic torque on B.
  // Note: Rotors 0 and 2 rotate one way and rotors 1 and 3 rotate the other.
  const auto tau_B = T::sym({Mx, My, Mz});

  // Calculate local celestial body's (Earth's) gravity force on B, expressed in
  // the Newtonian frame N (a.k.a the inertial or World frame).
  const auto Fgravity_N = T::sym({0, 0, -m_ * g_});

  // Extract roll-pitch-yaw angles (rpy) and their time-derivatives (rpyDt).
  const auto rpy = X_.rpy.positions();
  const auto rpyDt = X_.rpy.velocities();

  // Convert roll-pitch-yaw (rpy) orientation to the R_NB rotation matrix.
  const auto R_NB = rotation_matrix(rpy);

  // Calculate the net force on B, expressed in N.  Use Newton's law to
  // calculate a_NBcm_N (acceleration of B's center of mass, expressed in N).
  const auto Fnet_N = Fgravity_N + R_NB * Faero_B;
  const auto xyzDDt = Fnet_N / m_; // Equal to a_NBcm_N.

  // Use rpy and rpyDt to calculate B's angular velocity in N, expressed in B.
  const auto w_BN_B = angular_velocity_in_child_from_rpyDt(rpy, rpyDt);

  // To compute α (B's angular acceleration in N) due to the net moment 𝛕 on B,
  // rearrange Euler rigid body equation  𝛕 = I α + ω × (I ω)  and solve for α.
  const auto wIw = T::cross(w_BN_B, I_ * w_BN_B); // Expressed in B

  // Solve for α by solving linear equation I α = b
  // where b == 𝛕 - ω × (I ω)
  const auto b = tau_B - wIw;

  SX D, LT;
  std::vector<casadi_int> p;
  SX::ldl(I_, D, LT, p, false); // LDL factorization
  auto alpha_NB_B = SX::ldl_solve(b, D, LT, p); // Expressed in B

  const auto alpha_NB_N = R_NB * alpha_NB_B; // Expressed in N

  // Calculate the 2nd time-derivative of rpy.
  const auto rpyDDt = rpyDDT_from_rpyDt_and_angular_accel_in_parent(
      rpy, rpyDt, alpha_NB_N);

  return T::sym({X_.xyz.velocities(), X_.rpy.velocities(), xyzDDt, rpyDDt});
}

template <typename T>
PlanarQuadrotor<T>::PlanarQuadrotor()
    : PlanarQuadrotor(0.1, 0.1)
      {}

template <typename T>
PlanarQuadrotor<T>::PlanarQuadrotor(double m, double l_arg)
    {}

} // namespace quadrotor
