# pragma once

#include <cmath>
#include <stdexcept>

namespace rollpitchyaw {

inline constexpr double kGimbalLockToleranceCosPitchAngle = 0.008;

template <typename T>
T rotation_matrix(const T& rpy);

template <typename T>
T matrix_relating_angular_velocity_in_child_to_rpyDt(const T& rpy);

template <typename T>
T angular_velocity_in_child_from_rpyDt(const T& rpy, const T& rpyDt);

template <typename T>
T matrix_relating_rpyDt_to_angular_velocity_in_parent(const T& rpy);

template <typename T>
T matrix_relating_angular_velocity_in_parent_to_rpyDt(
    const T& rpy, const T& rpyDt);

template <typename T>
T rpyDDT_from_rpyDt_and_angular_accel_in_parent(
    const T& rpy, const T& rpyDt, const T& alpha_AD_A);

} // namespace rollpitchyaw
