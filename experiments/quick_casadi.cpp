#include <casadi/casadi.hpp>

using casadi::SX;
using casadi::DM;
using casadi::MX;
using casadi::Slice;
using std::cout;
using std::endl;

SX rotation_matrix(const SX& rpy) {
  auto& sin = SX::sin;
  auto& cos = SX::sin;
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
  const auto R_row1 = horzcat(Rxx, Rxy, Rxz);
  const auto R_row2 = horzcat(Ryx, Ryy, Ryz);
  const auto R_row3 = horzcat(Rzx, Rzy, Rzz);
  return vertcat(R_row1, R_row2, R_row3);
}

int main () {
  SX x = SX::sym("", 3);
  x(0) = 1;
  x(1) = 2;
  x(2) = 3;
  SX R = SX::eye(3);
  R(2,2) = 5;

  auto z = mtimes(R, x);

  cout << z.size() << endl;

  return 0;
}
