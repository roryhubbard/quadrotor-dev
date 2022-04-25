# pragma once

#include <cmath>
#include <stdexcept>
#include <casadi/casadi.hpp>

namespace rollpitchyaw {

using casadi::SX;

inline constexpr double kGimbalLockToleranceCosPitchAngle = 0.008;

SX rotation_matrix(const SX& rpy);

SX matrix_relating_angular_velocity_in_child_to_rpyDt(const SX& rpy);

SX angular_velocity_in_child_from_rpyDt(const SX& rpy, const SX& rpyDt);

SX matrix_relating_rpyDt_to_angular_velocity_in_parent(const SX& rpy);

SX matrix_relating_angular_velocity_in_parent_to_rpyDt(
    const SX& rpy, const SX& rpyDt);

SX rpyDDt_from_rpyDt_and_angular_accel_in_parent(
    const SX& rpy, const SX& rpyDt, const SX& alpha_AD_A);

} // namespace rollpitchyaw
