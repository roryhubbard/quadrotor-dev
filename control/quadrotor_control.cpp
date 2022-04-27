#include <iostream>
#include <cmath>
#include <cassert>
#include <casadi/casadi.hpp>
#include "quadrotor_plant.hpp"

namespace ca = casadi;
using quadrotor::Quadrotor;
using std::cout, std::endl;

template <typename T>
void append_vector(std::vector<T>& a, const std::vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
}

template <typename T>
std::vector<T> linspace(const T& a, const T& b, const std::size_t& N) {
  T h = (b - a) / static_cast<T>(N-1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
    *x = val;
  }
  return xs;
}

template<typename T>
std::vector<std::vector<T>> straight_line_trajectory_guess(
    const std::vector<T>& xi, const std::vector<T>& xf, const std::size_t& N) {

  const auto nstates = xi.size();

  std::vector<std::vector<T>> individual_traj;
  for (std::size_t i = 0; i < nstates; ++i) {
    const auto state_traj_i = linspace(xi[i], xf[i], N);
    individual_traj.push_back(state_traj_i);
  }

  std::vector<std::vector<T>> merged_traj(N);
  typename std::vector<std::vector<T>>::iterator x;
  std::size_t i = 0;
  for (x = merged_traj.begin(); x != merged_traj.end(); ++x) {
    std::vector<T> column(nstates);
    std::size_t j = 0;
    for (auto v = column.begin(); v != column.end(); ++v) {
      *v = std::move(individual_traj[j][i]);
      ++j;
    }
    *x = std::move(column);
    ++i;
  }

  assert(x==merged_traj.end());
  return merged_traj;
}

int main() {
  auto cf = Quadrotor();
  const auto X = cf.X();
  const auto U = cf.U();
  const auto dXdt = cf.ode();
  const auto nx = cf.nx();
  const auto nu = cf.nu();

  // objective
  ca::SX L = sumsqr(X);
  ca::SXDict dae = {{"x", X}, {"p", U}, {"ode", dXdt}, {"quad", L}};

  const double T = 1; // time horizon
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
  Xi[2] = -10;
  const std::vector<double> Xf(nx, 0.0);

  // bounds and guess for state
  const std::vector<double> lbX(nx, -20);
  const std::vector<double> ubX(nx,  20);
  const auto X_guess = straight_line_trajectory_guess(Xi, Xf, N+1);
  //const std::vector<double> X_guess(nx, 0.0);

  // bounds and guess for control
  const std::vector<double> lbU(nu, 0.0);
  const std::vector<double> ubU(nu, 4.0);
  const std::vector<double> U_guess(nu, 2.0);

  std::vector<ca::MX> w;            // decision variables
  std::vector<ca::MX> g;            // constraint function
  std::vector<double> lbw, ubw, w0; // decision variable bounds and guess
  std::vector<double> lbg, ubg;     // constraint bounds

  // start state
  auto Xk = ca::MX::sym("Xk_0", nx);
  w.push_back(Xk);
  append_vector(lbw, Xi);
  append_vector(ubw, Xi);
  append_vector(w0, X_guess[0]);
  //append_vector(w0, X_guess);

  ca::MX J = 0;

  for (auto k = 0; k < N; ++k) {
    const auto Uk = ca::MX::sym("Uk_" + std::to_string(k), nu);
    w.push_back(Uk);
    append_vector(lbw, lbU);
    append_vector(ubw, ubU);
    append_vector(w0, U_guess);

    // integrate till the end of the interval
    const auto Ik = F(ca::MXDict{{"x0", Xk}, {"p", Uk}});
    const auto Xk_end = Ik.at("xf");
    J += Ik.at("qf");

    const auto Xk = ca::MX::sym("Xk_" + std::to_string(k), nx);
    w.push_back(Xk);
    // continuity constraints
    g.push_back(Xk_end - Xk);

    if (k == N - 1) {
      append_vector(lbw, Xf);
      append_vector(ubw, Xf);
      append_vector(lbg, std::vector<double>(nx, -.1));
      append_vector(ubg, std::vector<double>(nx, 0.1));
    } else {
      append_vector(lbw, lbX);
      append_vector(ubw, ubX);
      append_vector(lbg, std::vector<double>(nx, 0.0));
      append_vector(ubg, std::vector<double>(nx, 0.0));
    }
    append_vector(w0, X_guess[k+1]);
    //append_vector(w0, X_guess);
  }

  // Total number of NLP variables
  const int NV = nx*(N+1) + nu*N;
  const auto w_vector = vertcat(w);
  assert(w_vector.numel()==NV);
  // NLP
  ca::MXDict nlp = {{"x", w_vector}, {"f", J}, {"g", vertcat(g)}};

  // Set options
  ca::Dict opts;
  opts["ipopt.tol"] = 1e-5;
  //opts["ipopt.max_iter"] = 100;

  // Create an NLP solver and buffers
  auto solver = nlpsol("nlpsol", "ipopt", nlp, opts);
  std::map<std::string, ca::DM> arg, res;

  // Bounds and initial guess
  arg["lbx"] = lbw;
  arg["ubx"] = ubw;
  arg["lbg"] = lbg;
  arg["ubg"] = ubg;
  arg["x0"] = w0;

  // Solve the problem
  res = solver(arg);

  // Optimal solution of the NLP
  std::vector<double> w_opt(res.at("x"));

  // Get the optimal state trajectory
  std::vector<double> x_opt(N+1), xdot_opt(N+1);
  std::vector<double> y_opt(N+1), ydot_opt(N+1);
  std::vector<double> z_opt(N+1), zdot_opt(N+1);
  std::vector<double> r_opt(N+1), rdot_opt(N+1);
  std::vector<double> p_opt(N+1), pdot_opt(N+1);
  std::vector<double> yaw_opt(N+1), yawdot_opt(N+1);

  for(int i=0; i<=N; ++i) {
    x_opt[i] = w_opt.at(        i*(nx+1));
    y_opt[i] = w_opt.at(   1  + i*(nx+1));
    z_opt[i] = w_opt.at(   2  + i*(nx+1));
    r_opt[i] = w_opt.at(   3  + i*(nx+1));
    p_opt[i] = w_opt.at(   4  + i*(nx+1));
    yaw_opt[i] = w_opt.at( 5  + i*(nx+1));
    xdot_opt[i] = w_opt.at(6  + i*(nx+1));
    ydot_opt[i] = w_opt.at(7  + i*(nx+1));
    zdot_opt[i] = w_opt.at(8  + i*(nx+1));
    rdot_opt[i] = w_opt.at(9  + i*(nx+1));
    pdot_opt[i] = w_opt.at(10 + i*(nx+1));
    yawdot_opt[i] = w_opt.at(11 + i*(nx+1));
  }

  cout << "x_opt = " << endl << x_opt << endl;
  cout << "y_opt = " << endl << y_opt << endl;
  cout << "z_opt = " << endl << z_opt << endl;
  cout << "r_opt = " << endl << r_opt << endl;
  cout << "p_opt = " << endl << p_opt << endl;
  cout << "yaw_opt = " << endl << yaw_opt << endl;

  // Get the optimal control
  std::vector<double> u1_opt(N);
  std::vector<double> u2_opt(N);
  std::vector<double> u3_opt(N);
  std::vector<double> u4_opt(N);

  for(int i=0; i<N; ++i) {
    u1_opt[i] = w_opt.at(nx +    i*(nx+1));
    u2_opt[i] = w_opt.at(nx + 1 +i*(nx+1));
    u3_opt[i] = w_opt.at(nx + 2 +i*(nx+1));
    u4_opt[i] = w_opt.at(nx + 3 +i*(nx+1));
  }

  cout << "u1_opt = " << endl << u1_opt << endl;
  cout << "u2_opt = " << endl << u2_opt << endl;
  cout << "u3_opt = " << endl << u3_opt << endl;
  cout << "u4_opt = " << endl << u4_opt << endl;

  return 0;
}
