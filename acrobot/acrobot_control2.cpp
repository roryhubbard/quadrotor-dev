#include <iostream>
#include <cmath>
#include <cassert>
#include <casadi/casadi.hpp>
#include <nlohmann/json.hpp>
#include "acrobot_plant.hpp"

using json = nlohmann::json;
namespace ca = casadi;
using acrobot::Acrobot;
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
std::vector<T> straight_line_trajectory_continuous(
    const std::vector<T>& xi, const std::vector<T>& xf, const std::size_t& N) {

  const auto nstates = xi.size();

  std::vector<std::vector<T>> individual_traj;
  for (std::size_t i = 0; i < nstates; ++i) {
    const auto state_traj_i = linspace(xi[i], xf[i], N);
    individual_traj.push_back(state_traj_i);
  }

  std::vector<std::vector<T>> merged_traj(N*nstates);
  auto x = merged_traj.begin();
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < nstates; ++j) {
      *x = std::move(individual_traj[j][i]);
    }
  }

  assert(x==merged_traj.end());
  return merged_traj;
}

template<typename T>
std::vector<std::vector<T>> straight_line_trajectory_stacked(
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

  return merged_traj;
}

ca::Function RK4_integrator(
    const double dt, const int N, const ca::SX& x,
    const ca::SX& u, const ca::SX& xdot, const ca::SX& L) {
  const int M = 1; // RK4 steps per interval
  const auto f = ca::Function("integrator", {x, u}, {xdot, L},
                                            {"x", "u"}, {"dx", "l"});

  auto X0 = ca::SX::sym("X0", 4);
  auto U = ca::SX::sym("U", 1);
  auto X = X0;
  auto Q = ca::SX::zeros(1);

  for (int j=0; j < M; ++j) {
    auto k1 = f(ca::SXDict{{"x", X}, {"u", U}});
    auto k2 = f(ca::SXDict{{"x", X + dt/2 * k1["dx"]}, {"u", U}});
    auto k3 = f(ca::SXDict{{"x", X + dt/2 * k2["dx"]}, {"u", U}});
    auto k4 = f(ca::SXDict{{"x", X + dt   * k3["dx"]}, {"u", U}});
    X += dt/6 * (k1["dx"] + 2*k2["dx"] + 2*k3["dx"] + k4["dx"]);
    Q += dt/6 * (k1["l"]  + 2*k2["l"]  + 2*k3["l"]  + k4["l"]);
  }

  return ca::Function("F", {X0, U}, {X, Q},
                           {"x0", "p"}, {"xf", "qf"});
}

ca::DMDict multiple_shooting(const Acrobot& acrobot,
                             const double& dt, const int& N) {
  const auto& X = acrobot.X();
  const auto& U = acrobot.U();
  const auto& ode = acrobot.ode();
  const auto nx = acrobot.nx();
  const auto nu = acrobot.nu();

  // Bounds and initial guess for the control
  const std::vector<double> u_min = {-10};
  const std::vector<double> u_max = {10};
  const std::vector<double> u_init = {0.0};

  // Bounds and initial guess for the state
  const std::vector<double> x0 = {   0, 0, 0, 0};
  const std::vector<double> xf = {M_PI, 0, 0, 0};
  const auto& x0_min = x0;
  const auto& x0_max = x0;
  const std::vector<double> x_min  = {-2*M_PI, -3*M_PI_4, -2*M_PI, -2*M_PI};
  const std::vector<double> x_max  = {2*M_PI,   3*M_PI_4,  2*M_PI,  2*M_PI};
  const auto& xf_min = xf;
  const auto& xf_max = xf;
  //const std::vector<std::vector<double>> x_init(N+1, {0, 0, 0, 0});
  const auto x_init = straight_line_trajectory_stacked(x0, xf, N+1);

  // objective
  ca::SX L = U*U;
  ca::SXDict dae = {{"x", X}, {"p", U}, {"ode", ode}, {"quad", L}};

  const auto F = RK4_integrator(dt, N, X, U, ode, L);

  // Total number of NLP variables
  const int NV = nx*(N+1) + nu*N;

  // Declare variable vector for the NLP
  auto V = ca::MX::sym("V", NV);

  // NLP variable bounds and initial guess
  std::vector<double> v_min, v_max, v_init;

  // Offset in V
  int offset=0;

  // State at each shooting node and control for each shooting interval
  std::vector<ca::MX> x, u;
  for(int k=0; k<N; ++k){
    // Local state
    x.push_back( V.nz(ca::Slice(offset, offset+nx)));
    if(k==0){
      v_min.insert(v_min.end(), x0_min.begin(), x0_min.end());
      v_max.insert(v_max.end(), x0_max.begin(), x0_max.end());
    } else {
      v_min.insert(v_min.end(), x_min.begin(), x_min.end());
      v_max.insert(v_max.end(), x_max.begin(), x_max.end());
    }
    v_init.insert(v_init.end(), x_init[k].begin(), x_init[k].end());
    offset += nx;

    // Local control
    u.push_back( V.nz(ca::Slice(offset,offset+nu)));
    v_min.insert(v_min.end(), u_min.begin(), u_min.end());
    v_max.insert(v_max.end(), u_max.begin(), u_max.end());
    v_init.insert(v_init.end(), u_init.begin(), u_init.end());
    offset += nu;
  }

  // State at end
  x.push_back(V.nz(ca::Slice(offset, offset+nx)));
  v_min.insert(v_min.end(), xf_min.begin(), xf_min.end());
  v_max.insert(v_max.end(), xf_max.begin(), xf_max.end());
  v_init.insert(v_init.end(), x_init.back().begin(), x_init.back().end());
  offset += nx;

  // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
  casadi_assert(offset==NV, "");

  // Objective function
  ca::MX J = 0;

  //Constraint function and bounds
  std::vector<ca::MX> g;

  // Loop over shooting nodes
  for(int k=0; k<N; ++k){
    // Create an evaluation node
    ca::MXDict I_out = F(ca::MXDict{{"x0", x[k]}, {"p", u[k]}});

    // Save continuity constraints
    g.push_back( I_out.at("xf") - x[k+1] );

    // Add objective function contribution
    J += I_out.at("qf");
  }

  // NLP
  ca::MXDict nlp = {{"x", V}, {"f", J}, {"g", vertcat(g)}};

  // Set options
  ca::Dict opts;
  opts["ipopt.tol"] = 1e-5;
  //opts["ipopt.max_iter"] = 100;

  // Create an NLP solver and buffers
  ca::Function solver = nlpsol("nlpsol", "ipopt", nlp, opts);
  std::map<std::string, ca::DM> arg, res;

  // Bounds and initial guess
  arg["lbx"] = v_min;
  arg["ubx"] = v_max;
  arg["lbg"] = 0;
  arg["ubg"] = 0;
  arg["x0"] = v_init;

  // Solve the problem
  return solver(arg);
}

int main() {
  auto acrobot = Acrobot();
  const auto nx = acrobot.nx();
  const auto nu = acrobot.nu();

  const double dt = .05; // control rate (Hz)
  const int N = 200; // # shooting nodes == # control ticks

  const auto res = multiple_shooting(acrobot, dt, N);

  // Optimal solution of the NLP
  std::vector<double> w_opt(res.at("x"));

  // Get the optimal state trajectory
  std::vector<double> q1_opt(N+1), q1dot_opt(N+1);
  std::vector<double> q2_opt(N+1), q2dot_opt(N+1);

  for(int i=0; i<=N; ++i) {
    q1_opt[i] = w_opt.at(i*(nx+1));
    q2_opt[i] = w_opt.at(1+i*(nx+1));
    q1dot_opt[i] = w_opt.at(2+i*(nx+1));
    q2dot_opt[i] = w_opt.at(3+i*(nx+1));
  }

  cout << "q1_opt = " << endl << q1_opt << endl;
  cout << "q2_opt = " << endl << q2_opt << endl;
  cout << "q1dot_opt = " << endl << q1dot_opt << endl;
  cout << "q2dot_opt = " << endl << q2dot_opt << endl;

  // Get the optimal control
  std::vector<double> u_opt(N);

  for(int i=0; i<N; ++i) {
    u_opt[i] = w_opt.at(nx+i*(nx+1));
  }

  cout << "u_opt = " << endl << u_opt << endl;

  json j;

  j["q1"] = q1_opt;
  j["q2"] = q2_opt;
  j["q1dot"] = q1dot_opt;
  j["q2dot"] = q2dot_opt;

  std::ofstream o("swingup_trajectory.json");
  o << std::setw(2) << j << std::endl;

  return 0;
}
