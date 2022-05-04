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

template<typename T>
T default_moment_of_inertia() {
  Slice all;
  T I = T::sym("I", 3, 3);
  I(0, all) = {0.0015, 0,      0};
  I(1, all) = {0, 0.0025,      0};
  I(2, all) = {0,      0, 0.0035};
  return I;
}

template<>
MX default_moment_of_inertia<MX>() {
  Slice all;
  MX I = MX::sym("I", 3, 3);
  I(0, 0) = 0.0015;
  I(0, 1) = 0;
  I(0, 2) = 0;
  I(1, 0) = 0;
  I(1, 1) = 0.0025;
  I(1, 2) = 0;
  I(2, 0) = 0;
  I(2, 1) = 0;
  I(2, 2) = 0.0035;
  return I;
}

int main () {
  MX V = MX::sym("V", 10);
  cout << V.nz(Slice(0, 3)) << endl;

  return 0;
}
