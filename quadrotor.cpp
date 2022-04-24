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
using std::vector;
using std::cout;
using std::endl;

template <typename T>
void append_vector(vector<T>& a, const vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
}

int main() {
  CrazyFlie cf = CrazyFlie();


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

  // bounds and guess for state
  const vector<double> lbX(nx, -inf);
  const vector<double> ubX(nx,  inf);
  const vector<double> X_guess(nx, 0.0);

  // bounds and guess for control
  const vector<double> lbU(nu, -1.0);
  const vector<double> ubU(nu, 1.0);
  const vector<double> U_guess(nu, 0.0);

  vector<MX> w;                // decision variables
  vector<MX> g;                // constraint function
  vector<double> lbw, ubw, w0; // decision variable bounds and guess

  MX Xk = MX::sym("Xk_0", nx);
  w.push_back(Xk);
  append_vector(lbw, Xi);
  append_vector(ubw, Xi);
  append_vector(w0, Xi);

  MX J = 0;

  for (int k=0; k<N; ++k) {
    MX Uk = MX::sym("Uk_" + std::to_string(k), nu);
    w.push_back(Uk);
    append_vector(lbw, lbU);
    append_vector(ubw, ubU);
    append_vector(w0, U_guess);

    // integrate till the end of the interval
    MXDict Ik = F(MXDict{{"x0", Xk}, {"p", Uk}});
    MX Xk_end = Ik.at("xf"); // can also use Ik["xf"]
    J += Ik.at("qf");

    MX Xk = MX::sym("Xk_" + std::to_string(k), nx);
    w.push_back(Xk);
    append_vector(lbw, lbX);
    append_vector(ubw, ubX);
    append_vector(w0, X_guess);

    // continuity constraints
    g.push_back(Xk_end - Xk);
  }

  // NLP
  MXDict nlp = {{"x", vertcat(w)}, {"f", J}, {"g", vertcat(g)}};

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
