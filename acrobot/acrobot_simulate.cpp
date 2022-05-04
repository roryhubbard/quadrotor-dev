#include <iostream>
#include <fstream>
#include <vector>
#include <casadi/casadi.hpp>
#include <nlohmann/json.hpp>
#include "acrobot_plant.hpp"

using json = nlohmann::json;
namespace ca = casadi;
using acrobot::Acrobot;

int main() {
  auto acrobot = Acrobot();
  const auto& X = acrobot.X();
  const auto& U = acrobot.U();
  const auto& dXdt = acrobot.ode();
  const auto nx = acrobot.nx();
  const auto nu = acrobot.nu();
  
  const double T = 2;
  const int N = 20;

  // objective
  ca::SX L = sumsqr(U);
  ca::SXDict dae = {{"x", X}, {"p", U}, {"ode", dXdt}, {"quad", L}};

  const bool use_cvodes = true;
  ca::Function F; // integrator

  if (use_cvodes) {
    // Create an integrator (CVodes)
    const double dt = T / N;
    F = integrator("integrator", "cvodes", dae, {{"t0", 0}, {"tf", dt}});
  } else {
    // Fixed step Runge-Kutta 4 integrator
    const int M = 1; // RK4 steps per interval
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

  std::vector<double> Xk(4, 0);
  Xk[0] = M_PI_4;
  const double Uk = 0;

  std::vector<std::vector<double>> X_simulated;
  X_simulated.push_back(Xk);

  for (int i=0; i<N; ++i) {
    // integrate till the end of the interval
    std::cout << Xk << std::endl;
    const auto Ik = F(ca::DMDict{{"x0", Xk}, {"p", Uk}});
    const std::vector<double> X_end(Ik.at("xf")); 
    X_simulated.push_back(X_end);
    Xk = X_end;
  }

  json j;

  j["X_simulated"] = X_simulated;

  std::ofstream o("pretty.json");
  o << std::setw(2) << j << std::endl;
}

