#pragma once

#include <casadi/casadi.hpp>

namespace systems {

using casadi::SX;

class System {

  public:
//    System();
//    virtual ~System();

    SX&                     X()       { return X_; }
    [[nodiscard]] const SX& X() const { return X_; }

    SX&                     U()       { return U_; }
    [[nodiscard]] const SX& U() const { return U_; }

    SX&                     ode()       { return ode_; }
    [[nodiscard]] const SX& ode() const { return ode_; }

    [[nodiscard]] int nx() const { return static_cast<int>(X_.numel()); }
    [[nodiscard]] int nu() const { return static_cast<int>(U_.numel()); }

  protected:
    virtual void set_ode() = 0;

  private:
    SX X_;   // state vector
    SX U_;   // control vector
    SX ode_; // dX/dt = f(x, u)
};

} // namespace systems

