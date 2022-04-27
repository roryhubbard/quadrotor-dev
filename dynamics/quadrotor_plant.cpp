#include "quadrotor_plant.hpp"

namespace quadrotor {

using casadi::SX;
namespace RPY = rollpitchyaw;

SX default_moment_of_inertia() {
  using casadi::Slice;
  Slice all;
  SX I = SX::sym("I", 3, 3);
  I(0, all) = {0.0015, 0,      0};
  I(1, all) = {0, 0.0025,      0};
  I(2, all) = {0,      0, 0.0035};
  return I;
}

Quadrotor::Quadrotor()
    : Quadrotor(0.775,  // m (kg)
                0.15,  // L (m)
                default_moment_of_inertia(),
                1.0,    // kF
                0.0245  // kM
                ) {}

Quadrotor::Quadrotor(double m, double L, SX I, double kF, double kM)
    : g_(9.81),
      m_(m),
      L_(L),
      I_(I),
      kF_(kF),
      kM_(kM),
      U_(SX::sym("U",4)),
      ode_(calculate_ode()) {}

SX Quadrotor::calculate_ode() {
  // https://andrew.gibiansky.com/downloads/pdf/Quadcopter%20Dynamics,%20Simulation,%20and%20Control.pdf
  // https://github.com/RobotLocomotion/drake/blob/master/examples/quadrotor/quadrotor_plant.cc#L62

  // For each rotor, calculate the Bz measure of its aerodynamic force on B.
  // Note: B is the quadrotor body and Bz is parallel to each rotor's spin axis.
  const auto uF_Bz = kF_ * U_;
  const auto thrust = sum1(uF_Bz);

  // Compute the net aerodynamic force on B (from the 4 rotors), expressed in B.
  auto Faero_B = SX::zeros(3);
  Faero_B(2) = thrust;

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
  auto tau_B = SX::sym("tau_B", 3);
  tau_B(0) = Mx, tau_B(1) = My, tau_B(2) = Mz;

  // Calculate local celestial body's (Earth's) gravity force on B, expressed in
  // the Newtonian frame N (a.k.a the inertial or World frame).
  auto Fgravity_N = SX::zeros(3);
  Fgravity_N(2) = -m_ * g_;

  // Extract roll-pitch-yaw angles (rpy) and their time-derivatives (rpyDt).
  const auto rpy = X_.rpy.positions();
  const auto rpyDt = X_.rpy.velocities();

  // Convert roll-pitch-yaw (rpy) orientation to the R_NB rotation matrix.
  const auto R_NB = RPY::rotation_matrix(rpy);

  // Calculate the net force on B, expressed in N.  Use Newton's law to
  // calculate a_NBcm_N (acceleration of B's center of mass, expressed in N).
  const auto Fnet_N = Fgravity_N + mtimes(R_NB, Faero_B);
  const auto xyzDDt = Fnet_N / m_; // Equal to a_NBcm_N.

  // Use rpy and rpyDt to calculate B's angular velocity in N, expressed in B.
  const auto w_BN_B = RPY::angular_velocity_in_child_from_rpyDt(rpy, rpyDt);

  // To compute α (B's angular acceleration in N) due to the net moment 𝛕 on B,
  // rearrange Euler rigid body equation  𝛕 = I α + ω × (I ω)  and solve for α.
  const auto wIw = cross(w_BN_B, mtimes(I_, w_BN_B)); // Expressed in B

  // Solve for α by solving linear equation I α = b
  // where b == 𝛕 - ω × (I ω)
  const auto b = tau_B - wIw;

  SX D, LT;
  std::vector<casadi_int> p;
  ldl(I_, D, LT, p, false); // LDL factorization
  const auto alpha_NB_B = ldl_solve(b, D, LT, p); // Expressed in B

  const auto alpha_NB_N = mtimes(R_NB, alpha_NB_B); // Expressed in N

  // Calculate the 2nd time-derivative of rpy.
  const auto rpyDDt = RPY::rpyDDt_from_rpyDt_and_angular_accel_in_parent(
      rpy, rpyDt, alpha_NB_N);

  return vertcat(X_.velocities(), xyzDDt, rpyDDt);
}

PlanarQuadrotor::PlanarQuadrotor()
    : PlanarQuadrotor(0.1, 0.1)
      {}

PlanarQuadrotor::PlanarQuadrotor(double m, double l_arg)
    {}

} // namespace quadrotor

