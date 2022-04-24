#include <casadi/casadi.hpp>

using casadi::SX;
using casadi::MX;
using casadi::Linsol;
using std::cout;
using std::endl;

template <typename T>
void default_moment_of_inertia(T& I) {
  using casadi::Slice;
  Slice all;
  I(0, all) = {2.0, 0, 0};
  I(1, all) = {0, 1.0, 0};
  I(2, all) = {0, 0, 1.0};
}

template <>
void default_moment_of_inertia(MX& I) {
  I(0, 0) = 2.0;
  I(1, 1) = 1.0;
  I(2, 2) = 1.0;
}

int main () {
  auto A = SX::sym("", 3, 3);
  auto& c = SX::cos;

  cout << c(A) << endl;
  cout << 1 / 2.1 << endl;
 // default_moment_of_inertia(A);

 // auto b = SX::sym("", 3, 1);
//  auto b = SX::ones(3);

//  Linsol F("F", "ma27", A.sparsity());
//  MX r = F.solve(A, b);

//  SX D, LT;
//  std::vector<casadi_int> p;
//  SX::ldl(A, D, LT, p, false);
//  SX r = SX::ldl_solve(b, D, LT, p);
//
//  cout << A << endl;
//  cout << D << endl;
//  cout << LT << endl;
//  cout << r << endl;
//  cout << r(0) << endl;
//  cout << r(1) << endl;
//  cout << r(2) << endl;


  return 0;
}
