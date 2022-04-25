#include <iostream>
#include <cmath>
#include <casadi/casadi.hpp>
#include "quadrotor_plant.hpp"

namespace ca = casadi;
using quadrotor::QuadRotor;
using std::cout, std::endl;

template <typename T>
void append_vector(std::vector<T>& a, const std::vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
}

int main() {
  auto cf = QuadRotor();
  const auto X = cf.X();
  const auto U = cf.U();
  const auto dXdt = cf.ode();
  const auto nx = cf.nx();
  const auto nu = cf.nu();

  // objective
  ca::SX L = sumsqr(X);
  ca::SXDict dae = {{"x", X}, {"p", U}, {"ode", dXdt}, {"quad", L}};

  const double T = 4; // time horizon
  const int N = 10; // # shooting nodes == # control ticks

  const bool use_cvodes = true;
  ca::Function F; // integrator

  if (use_cvodes) {
    // Create an integrator (CVodes)
    F = integrator("integrator", "cvodes", dae, {{"t0", 0}, {"tf", T / N}});
  } else {
    // Fixed step Runge-Kutta 4 integrator
    const int M = 4; // RK4 steps per interval
    const double DT = T / N / M;
    const auto f = ca::Function("integrator", {X, U}, {dXdt, L},
                                              {"x", "u"}, {"dx", "l"});

    auto X0 = ca::SX::sym("X0", nx);
    auto U = ca::SX::sym("U", nu);
    auto X = X0;
    auto Q = ca::SX::sym("Q");

    for (int j=0; j < M; ++j) {
      auto k1 = f(ca::SXDict{{"x", X}, {"u", U}});
      auto k2 = f(ca::SXDict{{"x", X + DT/2 * k1["dx"]}, {"u", U}});
      auto k3 = f(ca::SXDict{{"x", X + DT/2 * k2["dx"]}, {"u", U}});
      auto k4 = f(ca::SXDict{{"x", X + DT   * k3["dx"]}, {"u", U}});
      X += DT/6 * (k1["dx"] + 2*k2["dx"] + 2*k3["dx"] + k4["dx"]);
      Q += DT/6 * (k1["l"]  + 2*k2["l"]  + 2*k3["l"]  + k4["l"]);
    }

    F = ca::Function("F", {X0, U}, {X, Q},
                          {"x0", "p"}, {"xf", "qf"});
  }

  // initial and final state
  std::vector<double> Xi(nx, 0.0);
  Xi[0] = -1;

  // bounds and guess for state
  const std::vector<double> lbX(nx, -ca::inf);
  const std::vector<double> ubX(nx,  ca::inf);
  const std::vector<double> X_guess(nx, 0.0);

  // bounds and guess for control
  const std::vector<double> lbU(nu, 0.0);
  const std::vector<double> ubU(nu, 1.0);
  const std::vector<double> U_guess(nu, 0.0);

  std::vector<ca::MX> w;            // decision variables
  std::vector<ca::MX> g;            // constraint function
  std::vector<double> lbw, ubw, w0; // decision variable bounds and guess

  auto Xk = ca::MX::sym("Xk_0", nx);
  w.push_back(Xk);
  append_vector(lbw, Xi);
  append_vector(ubw, Xi);
  append_vector(w0, Xi);

  ca::MX J = 0;

  for (int k=0; k<N; ++k) {
    auto Uk = ca::MX::sym("Uk_" + std::to_string(k), nu);
    w.push_back(Uk);
    append_vector(lbw, lbU);
    append_vector(ubw, ubU);
    append_vector(w0, U_guess);

    // integrate till the end of the interval
    auto Ik = F(ca::MXDict{{"x0", Xk}, {"p", Uk}});
    auto Xk_end = Ik.at("xf");
    J += Ik.at("qf");

    auto Xk = ca::MX::sym("Xk_" + std::to_string(k), nx);
    w.push_back(Xk);
    append_vector(lbw, lbX);
    append_vector(ubw, ubX);
    append_vector(w0, X_guess);

    // continuity constraints
    g.push_back(Xk_end - Xk);
  }

   // NLP
   ca::MXDict nlp = {{"x", vertcat(w)}, {"f", J}, {"g", vertcat(g)}};

   // Set options
   ca::Dict opts;
//   opts["ipopt.tol"] = 1e-5;
//   opts["ipopt.max_iter"] = 100;

   // Create an NLP solver and buffers
   auto solver = nlpsol("nlpsol", "ipopt", nlp, opts);
   std::map<std::string, ca::DM> arg, res;

   // Bounds and initial guess
   arg["lbx"] = lbw;
   arg["ubx"] = ubw;
   arg["lbg"] = 0.0;
   arg["ubg"] = 0.0;
   arg["x0"] = w0;

   // Solve the problem
   res = solver(arg);

   // Optimal solution of the NLP
   std::vector<double> w_opt(res.at("x"));

   // Get the optimal state trajectory
   std::vector<double> y_opt(N+1), ydot_opt(N+1);
   std::vector<double> z_opt(N+1), zdot_opt(N+1);
   std::vector<double> theta_opt(N+1), thetadot_opt(N+1);

   for(int i=0; i<=N; ++i) {
     y_opt[i] = w_opt.at(i*(nx+1));
     ydot_opt[i] = w_opt.at(1+i*(nx+1));
     z_opt[i] = w_opt.at(2+i*(nx+1));
     zdot_opt[i] = w_opt.at(3+i*(nx+1));
     theta_opt[i] = w_opt.at(4+i*(nx+1));
     thetadot_opt[i] = w_opt.at(5+i*(nx+1));
   }

   cout << "y_opt = " << endl << y_opt << endl;
   cout << "theta_opt = " << endl << theta_opt << endl;

   // Get the optimal control
   std::vector<double> ur_opt(N);
   std::vector<double> ul_opt(N);

   for(int i=0; i<N; ++i) {
     ur_opt[i] = w_opt.at(nx+i*(nx+1));
     ul_opt[i] = w_opt.at(nx+1+i*(nx+1));
   }

   cout << "ur_opt = " << endl << ur_opt << endl;
   cout << "ul_opt = " << endl << ur_opt << endl;

  return 0;
}
