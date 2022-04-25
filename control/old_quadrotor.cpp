#include <iostream>
#include <cmath>
#include <casadi/casadi.hpp>

using casadi::SX;
using casadi::SXDict;
using casadi::MX;
using casadi::MXDict;
using casadi::DM;
using casadi::Dict;
using casadi::inf;
using casadi::Function;
using casadi::Slice;
using std::vector;
using std::cout;
using std::endl;

struct CrazyFlie {
  // constant parameters, crazyflie 2.1 specs
  const double m; // mass (kg)
  const double l; // length (m)
  const double I; // rotational inertia, uniform bar
  const double g = 9.8; // gravity

  explicit CrazyFlie(double m=0.027, double l=0.092):
    m(m),
    l(l),
    I(1./12. * m * std::pow(l, 2))
    {}
};

template <typename T>
void append_vector(vector<T>& a, const vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
}

int main() {
  // planar quadrotor first order individual states
  SX y = SX::sym("y");
  SX ydot = SX::sym("ydot");
  SX z = SX::sym("z");
  SX zdot = SX::sym("zdot");
  SX theta = SX::sym("theta");
  SX thetadot = SX::sym("thetadot");

  // thrust of propellors
  SX ur = SX::sym("ur"); // right
  SX ul = SX::sym("ul"); // left

  // state and control vectors
  SX X = vertcat(y, ydot, z, zdot, theta, thetadot);
  SX U = vertcat(ur, ul);

  const int nx = static_cast<int>(X.size1());
  const int nu = static_cast<int>(U.size1());

  CrazyFlie cf = CrazyFlie();
  // ODE: dX/dt = f(x, u)
  SX dXdt = vertcat(
      ydot,
      (-(ur + ul) * SX::sin(theta)) / cf.m,
      zdot,
      ((ur + ul) * SX::cos(theta)) / cf.m - cf.g,
      thetadot,
      cf.l * (ur - ul) / (2 * cf.I)
  );

  // objective
  SX L = SX::sumsqr(X);
  SXDict dae = {{"x", X}, {"p", U}, {"ode", dXdt}, {"quad", L}};

  const double T = 10; // time horizon
  const int N = 10; // # shooting nodes == # control ticks

  const bool use_cvodes = true;
  Function F; // integrator

  if (use_cvodes) {
    // Create an integrator (CVodes)
    F = integrator("integrator", "cvodes", dae, {{"t0", 0}, {"tf", T / N}});
  } else {
    // Fixed step Runge-Kutta 4 integrator
    const int M = 4; // RK4 steps per interval
    const double DT = T / N / M;
    Function f = Function("integrator", {X, U}, {dXdt, L},
                                        {"x", "u"}, {"dx", "l"});

    SX X0 = SX::sym("X0", nx);
    SX U = SX::sym("U", nu);
    SX X = X0;
    SX Q = SX::sym("Q");

    for (int j=0; j < M; ++j) {
      SXDict k1 = f(SXDict{{"x", X}, {"u", U}});
      SXDict k2 = f(SXDict{{"x", X + DT/2 * k1["dx"]}, {"u", U}});
      SXDict k3 = f(SXDict{{"x", X + DT/2 * k2["dx"]}, {"u", U}});
      SXDict k4 = f(SXDict{{"x", X + DT   * k3["dx"]}, {"u", U}});
      X += DT/6 * (k1["dx"] + 2*k2["dx"] + 2*k3["dx"] + k4["dx"]);
      Q += DT/6 * (k1["l"]  + 2*k2["l"]  + 2*k3["l"]  + k4["l"]);
    }

    F = Function("F", {X0, U}, {X, Q},
                      {"x0", "p"}, {"xf", "qf"});
  }

  // initial and final state
  const vector<double> Xi = {-10, 0, 0, 0, 0, 0};
  const vector<double> Xf = { 0, 0, 0, 0, 0, 0};

  // bounds and guess for state
  const vector<double> lbX(nx, -inf);
  const vector<double> ubX(nx,  inf);
  const vector<double> X_guess(nx, 0.0);

  // bounds and guess for control
  const vector<double> lbU(nu, -1.0);
  const vector<double> ubU(nu, 1.0);
  const vector<double> U_guess(nu, 0.0);

  // Total number of NLP variables
  int nw = nx*(N+1) + nu*N;

  // Declare variable vector for the NLP
  MX w = MX::sym("w", nw);

  // NLP variable bounds and initial guess
  vector<double> w0, lbw, ubw;

  // Offset in V
  int offset=0;

  // State at each shooting node and control for each shooting interval
  vector<MX> Xk, Uk;
  for (int k=0; k<N; ++k) {
    // Local state
    Xk.push_back(w.nz(Slice(offset, offset+nx)));
    if (k == 0) {
      lbw.insert(lbw.end(), Xi.begin(), Xi.end());
      ubw.insert(ubw.end(), Xi.begin(), Xi.end());
    } else {
      lbw.insert(lbw.end(), lbX.begin(), lbX.end());
      ubw.insert(ubw.end(), ubX.begin(), ubX.end());
    }
    w0.insert(w0.end(), X_guess.begin(), X_guess.end());
    offset += nx;

    // Local control
    Uk.push_back(w.nz(Slice(offset, offset+nu)));
    lbw.insert(lbw.end(), lbU.begin(), lbU.end());
    ubw.insert(ubw.end(), ubU.begin(), ubU.end());
    w0.insert(w0.end(), U_guess.begin(), U_guess.end());
    offset += nu;
  }

  // State at end
  Xk.push_back(w.nz(Slice(offset, offset+nx)));
  lbw.insert(lbw.end(), Xf.begin(), Xf.end());
  ubw.insert(ubw.end(), Xf.begin(), Xf.end());
  w0.insert(w0.end(), X_guess.begin(), X_guess.end());
  offset += nx;

  casadi_assert(offset==nw, "");

  // Objective function
  MX J = 0;

  //Constraint function and bounds
  vector<MX> g;

  // Loop over shooting nodes
  for (int k=0; k<N; ++k) {
    // Create an evaluation node
    MXDict I_out = F(MXDict{{"x0", Xk[k]}, {"p", Uk[k]}});

    // Save continuity constraints
    g.push_back( I_out.at("xf") - Xk[k+1] );

    // Add objective function contribution
    J += I_out.at("qf");
  }

  // NLP
  MXDict nlp = {{"x", w}, {"f", J}, {"g", vertcat(g)}};

  // Set options
  Dict opts;
  opts["ipopt.tol"] = 1e-5;
  opts["ipopt.max_iter"] = 100;

  // Create an NLP solver and buffers
  Function solver = nlpsol("nlpsol", "ipopt", nlp, opts);
  std::map<std::string, DM> arg, res;

  // Bounds and initial guess
  arg["lbx"] = lbw;
  arg["ubx"] = ubw;
  arg["lbg"] = 0.0;
  arg["ubg"] = 0.0;
  arg["x0"] = w0;

  // Solve the problem
  res = solver(arg);

  // Optimal solution of the NLP
  vector<double> w_opt(res.at("x"));

  // Get the optimal state trajectory
  vector<double> y_opt(N+1), ydot_opt(N+1);
  vector<double> z_opt(N+1), zdot_opt(N+1);
  vector<double> theta_opt(N+1), thetadot_opt(N+1);

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
  vector<double> ur_opt(N);
  vector<double> ul_opt(N);

  for(int i=0; i<N; ++i) {
    ur_opt[i] = w_opt.at(nx+i*(nx+1));
    ul_opt[i] = w_opt.at(nx+1+i*(nx+1));
  }

  cout << "ur_opt = " << endl << ur_opt << endl;
  cout << "ul_opt = " << endl << ur_opt << endl;

  return 0;
}
