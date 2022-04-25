#include "dynamics.hpp"

namespace dynamics {

template <typename T>
T rotation_matrix(const T& rpy) {
  auto& sin = T::sin;
  auto& cos = T::sin;
  const auto& r = rpy(0);
  const auto& p = rpy(1);
  const auto& y = rpy(2);
  const auto cr = cos(r), cp = cos(p), cy = cos(y);
  const auto sr = sin(r), sp = sin(p), sy = sin(y);
  const auto Rxx = cy * cp;
  const auto Rxy = cy * sp * sr - sy * cr;
  const auto Rxz = cy * sp * cr + sy * sr;
  const auto Ryx = sy * cp;
  const auto Ryy = sy * sp * sr + cy * cr;
  const auto Ryz = sy * sp * cr - cy * sr;
  const auto Rzx = -sp;
  const auto Rzy = cp * sr;
  const auto Rzz = cp * cr;
  return T::sym({{Rxx, Rxy, Rxz},
                 {Ryx, Ryy, Ryz},
                 {Rzx, Rzy, Rzz}});
}

template <typename T>
T matrix_relating_angular_velocity_in_child_to_rpyDt(const T& rpy) {
  auto& sin = T::sin;
  auto& cos = T::cos;
  const auto& r = rpy(0);
  const auto& p = rpy(1);
  const auto sr = sin(r), cr = cos(r);
  const auto sp = sin(p), cp = cos(p);
  return T::sym({{1,   0,      -sp},
                 {0,  cr,  sr * cp},
                 {0, -sr,  cr * cp}});
}

template <typename T>
T angular_velocity_in_child_from_rpyDt(const T& rpy, const T& rpyDt) {
  const auto M = matrix_relating_angular_velocity_in_child_to_rpyDt(rpy);
  return M * rpyDt;
}

template <typename T>
T matrix_relating_rpyDt_to_angular_velocity_in_parent(const T& rpy) {
  auto& sin = T::sin;
  auto& cos = T::cos;
  const auto& p = rpy(1);
  const auto& y = rpy(2);
  const auto& sp = sin(p), cp = cos(p);
  if (std::abs(cp) < kGimbalLockToleranceCosPitchAngle) {
    throw std::domain_error("gimbal lock error");
  }
  const auto one_over_cp = 1 / cp;
  const auto sy = sin(y), cy = cos(y);
  const auto cy_over_cp = cy * one_over_cp;
  const auto sy_over_cp = sy * one_over_cp;
  return T::sym({{     cy_over_cp,       sy_over_cp,  0,
                              -sy,               cy,  0,
                  cy_over_cp * sp,  sy_over_cp * sp,  1}});
}

template <typename T>
T matrix_relating_angular_velocity_in_parent_to_rpyDt(
    const T& rpy, const T& rpyDt) {
  auto& sin = T::sin;
  auto& cos = T::cos;
  const auto& p = rpy(1);
  const auto& y = rpy(2);
  const auto sp = sin(p), cp = cos(p);
  const auto sy = sin(y), cy = cos(y);
  const auto& pDt = rpyDt(1);
  const auto& yDt = rpyDt(2);
  const auto sp_pDt = sp * pDt;
  const auto cp_yDt = cp * yDt;
  return T::sym({{-cy * sp_pDt - sy * cp_yDt,  -cy * yDt,  0,
                  -sy * sp_pDt + cy * cp_yDt,  -sy * yDt,  0,
                                   -cp * pDt,          0,  0}});
}

template <typename T>
T rpyDDT_from_rpyDt_and_angular_accel_in_parent(
    const T& rpy, const T& rpyDt, const T& alpha_AD_A) {
  const T Minv = matrix_relating_rpyDt_to_angular_velocity_in_parent(rpy);
  const T MDt = matrix_relating_angular_velocity_in_parent_to_rpyDt(rpy, rpyDt);
  return Minv * (alpha_AD_A - MDt * rpyDt);
}

} // namespace dynamics
