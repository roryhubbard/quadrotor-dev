#include "roll_pitch_yaw.hpp"

namespace rollpitchyaw {

SX rotation_matrix(const SX& rpy) {
  auto& sin = SX::sin;
  auto& cos = SX::sin;
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
  const auto R_row1 = horzcat(Rxx, Rxy, Rxz);
  const auto R_row2 = horzcat(Ryx, Ryy, Ryz);
  const auto R_row3 = horzcat(Rzx, Rzy, Rzz);
  return vertcat(R_row1, R_row2, R_row3);
}

SX matrix_relating_angular_velocity_in_child_to_rpyDt(const SX& rpy) {
  auto& sin = SX::sin;
  auto& cos = SX::cos;
  const auto& r = rpy(0);
  const auto& p = rpy(1);
  const auto sr = sin(r), cr = cos(r);
  const auto sp = sin(p), cp = cos(p);
  const auto M_row1 = horzcat(1,   0,      -sp);
  const auto M_row2 = horzcat(0,  cr,  sr * cp);
  const auto M_row3 = horzcat(0, -sr,  cr * cp);
  return vertcat(M_row1, M_row2, M_row3);
}

SX angular_velocity_in_child_from_rpyDt(const SX& rpy, const SX& rpyDt) {
  const auto M = matrix_relating_angular_velocity_in_child_to_rpyDt(rpy);
  return mtimes(M, rpyDt);
}

SX matrix_relating_rpyDt_to_angular_velocity_in_parent(const SX& rpy) {
  auto& sin = SX::sin;
  auto& cos = SX::cos;
  const auto& p = rpy(1);
  const auto& y = rpy(2);
  const auto& sp = sin(p), cp = cos(p);
  const auto one_over_cp = 1 / cp;
  const auto sy = sin(y), cy = cos(y);
  const auto cy_over_cp = cy * one_over_cp;
  const auto sy_over_cp = sy * one_over_cp;
  const auto M_row1 = horzcat(     cy_over_cp,       sy_over_cp,  0);
  const auto M_row2 = horzcat(            -sy,               cy,  0);
  const auto M_row3 = horzcat(cy_over_cp * sp,  sy_over_cp * sp,  1);
  return vertcat(M_row1, M_row2, M_row3);
}

SX matrix_relating_angular_velocity_in_parent_to_rpyDt(
    const SX& rpy, const SX& rpyDt) {
  auto& sin = SX::sin;
  auto& cos = SX::cos;
  const auto& p = rpy(1);
  const auto& y = rpy(2);
  const auto sp = sin(p), cp = cos(p);
  const auto sy = sin(y), cy = cos(y);
  const auto& pDt = rpyDt(1);
  const auto& yDt = rpyDt(2);
  const auto sp_pDt = sp * pDt;
  const auto cp_yDt = cp * yDt;
  const auto M_row1 = horzcat(-cy * sp_pDt - sy * cp_yDt,  -cy * yDt,  0);
  const auto M_row2 = horzcat(-sy * sp_pDt + cy * cp_yDt,  -sy * yDt,  0);
  const auto M_row3 = horzcat(                 -cp * pDt,          0,  0);
  return vertcat(M_row1, M_row2, M_row3);
}

SX rpyDDt_from_rpyDt_and_angular_accel_in_parent(
    const SX& rpy, const SX& rpyDt, const SX& alpha_AD_A) {
  const SX Minv = matrix_relating_rpyDt_to_angular_velocity_in_parent(rpy);
  const SX MDt = matrix_relating_angular_velocity_in_parent_to_rpyDt(rpy, rpyDt);
  return mtimes(Minv,alpha_AD_A - mtimes(MDt, rpyDt));
}

} // namespace rollpitchyaw
