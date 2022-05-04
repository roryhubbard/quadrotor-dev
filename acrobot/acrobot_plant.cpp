#include "acrobot_plant.hpp"

namespace acrobot {

using casadi::SX;
using casadi::Slice;

Acrobot::Acrobot()
    : g_(9.81),
      m1_(1.0),
      m2_(1.0),
      l1_(1.0),
      l2_(1.0),
      I1_(1./3.),
      I2_(1./3.) {
    this->X() = SX::sym("X", 4);
    this->U() = SX::sym("U", 1);
    set_ode();
  }


void Acrobot::set_ode() {
  // http://underactuated.mit.edu/acrobot.html
  const auto& X = this->X();
  const auto l1c = l1_ / 2;
  const auto l2c = l2_ / 2;
  const auto ml1l2c = m2_ * l1_ * (l2c/2);

  auto M = SX::sym("M", 2, 2);
  M(0, 0) = I1_ + I2_ + m2_ * std::pow(l1_, 2) + 2 * ml1l2c * cos(X(1));
  M(0, 1) = I2_ + ml1l2c * cos(X(1));
  M(1, 0) = I2_ + ml1l2c * cos(X(1));
  M(1, 1) = I2_;

  auto C = SX::sym("C", 2, 2);
  C(0, 0) = -2 * ml1l2c * sin(X(1)) * X(3);
  C(0, 1) =     -ml1l2c * sin(X(1)) * X(3);
  C(1, 0) =      ml1l2c * sin(X(1)) * X(2);
  C(1, 1) = 0;

  auto tau = SX::sym("tau", 2, 1);
  tau(0) = -m1_*g_*l1c*sin(X(0)) - m2_*g_*(l1_*sin(X(0)) + l2c*sin(X(0) + X(1)));
  tau(1) = -m2_*g_*l2c*sin(X(0) + X(1));

  auto B = SX::sym("B", 2, 1);
  B(0) = 0;
  B(1) = 1;

  // solve for qddot: M*qddot + C*qdot = tau + B*U
  const auto b = tau + B * this->U() - mtimes(C, X(Slice(2, 4)));
  //SX D, LT;
  //std::vector<casadi_int> p;
  //ldl(M, D, LT, p, false); // LDL factorization
  //const auto qddot = ldl_solve(b, D, LT, p);
  const auto qddot = mtimes(inv(M), b);

  this->ode() = vertcat(X(Slice(2, 4)), qddot);
}

} // namespace acrobot
