#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <casadi/casadi.hpp>

using Eigen::Matrix3d;
using Eigen::Vector3d;
using Eigen::Ref;
using casadi::SX;
using casadi::DM;
using casadi::MX;
using std::cout;
using std::endl;

void default_moment_of_inertia(SX& I) {
  using casadi::Slice;
  Slice all;
  I(0, all) = {2.0, 0, 0};
  I(1, all) = {0, 1.0, 0};
  I(2, all) = {0, 0, 1.0};
}

void default_moment_of_inertia(Ref<Matrix3d> I) {
  I << 2.0,   0,  0.,
         0, 1.0,   0,
         0,   0, 1.0;
}

void default_moment_of_inertia(MX& I) {
  I(0, 0) = 2.0;
  I(1, 1) = 1.0;
  I(2, 2) = 1.0;
}

int main () {
  SX A_casadi = SX::sym("", 3, 3);
  default_moment_of_inertia(A_casadi);

  Matrix3d A_eigen;
  default_moment_of_inertia(A_eigen);

  const SX b_casadi = SX::ones(3);
  Vector3d b_eigen = Vector3d::Ones();

//  Linsol F("F", "ma27", A.sparsity());
//  MX r = F.solve(A, b);

  SX D, LT;
  std::vector<casadi_int> p;
  SX::ldl(A_casadi, D, LT, p, false);
  const SX x_casadi = SX::ldl_solve(b_casadi, D, LT, p);

  const Vector3d x_eigen = A_eigen.ldlt().solve(b_eigen);

  cout << x_casadi << endl;
  cout << x_eigen << endl;

  return 0;
}
