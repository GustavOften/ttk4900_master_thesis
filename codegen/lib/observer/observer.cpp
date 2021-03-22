//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  observer.cpp
//
//  Code generation for function 'observer'
//


// Include files
#include "observer.h"
#include "mod.h"
#include "observer_data.h"
#include "observer_initialize.h"
#include "observer_rtwutil.h"
#include "ppval.h"
#include "rt_nonfinite.h"
#include "test_controller.h"
#include <cmath>
#include <cstring>

// Variable Definitions
static struct_T spline_phi;
static bool spline_phi_not_empty;
static double q_hat[2];
static double dq_hat[2];
static double t_prev;

// Function Definitions
void observer(double t, double theta, double varphi, double u, double q[2],
              double dq[2])
{
  double phi;
  double dd_s;
  double d_sf;
  double d_gamma;
  double delta;
  double d_delta;
  double dd_delta;
  double scale;
  double a;
  double b_t;
  double ddd_delta_v[3];
  double m_12;
  double norm_d_delta_v;
  double delta_v[3];
  double d;
  double d1;
  double d_delta_v[3];
  double dd_delta_v[3];
  double tau[3];
  int r1;
  double d_tau_tmp[9];
  double b_dd_delta_v[9];
  double d_e_phi[3];
  double e_phi[3];
  double d_tau[3];
  double dd_tau[3];
  double b_gamma;
  static const signed char b[3] = { 1, 0, 0 };

  static const signed char b_b[3] = { 0, 1, 0 };

  double d_s_tmp;
  double M[4];
  static const signed char b_iv[3] = { 0, 0, 1 };

  double B[2];
  int r2;
  if (!isInitialized_observer) {
    observer_initialize();
  }

  if (!spline_phi_not_empty) {
    q_hat[0] = theta;
    q_hat[1] = varphi;
    std::memcpy(&spline_phi.breaks[0], &dv[0], 500U * sizeof(double));
    std::memcpy(&spline_phi.coefs[0], &dv1[0], 1996U * sizeof(double));
    spline_phi_not_empty = true;
  }

  // n_varphi = solutions.n_varphi;
  // Constant for describing the butterfly frame
  // Constant for describing the butterfly frame
  // Gravity
  // 5.48e-7;%5.48e-7; %Mass moment of inertia of the ball
  // Mass moment of inertia of the frame, from article.
  // Mass of the ball
  // R_b = sqrt(16.55e-3^2-12.5e-3^2); From paper.
  // Measured on the robot
  phi = ppval(spline_phi.breaks, spline_phi.coefs, b_mod(q_hat[1]));

  // Need to do this becaouse the simulink coder is strange...
  //     %% e_phi
  dd_s = std::sin(phi);
  d_sf = std::cos(phi);

  //     %% Delta
  d_gamma = std::cos(2.0 * phi);
  delta = 0.1095 - 0.0405 * d_gamma;
  phi = std::sin(2.0 * phi);
  d_delta = 0.081 * phi;
  dd_delta = 0.162 * d_gamma;
  scale = -0.324 * phi;
  a = 2.0 * d_delta;
  phi = 3.0 * dd_delta;
  d_gamma = 3.0 * d_delta;
  b_t = delta * -d_sf;
  ddd_delta_v[0] = ((scale * dd_s + phi * d_sf) + d_gamma * -dd_s) + b_t;
  m_12 = delta * dd_s;
  ddd_delta_v[1] = ((scale * d_sf + phi * -dd_s) + d_gamma * -d_sf) + m_12;
  ddd_delta_v[2] = ((scale * 0.0 + phi * 0.0) + d_gamma * 0.0) + delta * 0.0;
  norm_d_delta_v = std::sqrt(d_delta * d_delta + delta * delta);

  //     %% Tau
  delta_v[0] = m_12;
  d = delta * d_sf;
  d1 = d_delta * dd_s + d;
  d_delta_v[0] = d1;
  dd_delta_v[0] = (dd_delta * dd_s + a * d_sf) + delta * -dd_s;
  tau[0] = d1 / norm_d_delta_v;
  delta_v[1] = d;
  d1 = d_delta * d_sf + delta * -dd_s;
  d_delta_v[1] = d1;
  dd_delta_v[1] = (dd_delta * d_sf + a * -dd_s) + b_t;
  tau[1] = d1 / norm_d_delta_v;
  delta_v[2] = delta * 0.0;
  d1 = d_delta * 0.0 + delta * 0.0;
  d_delta_v[2] = d1;
  dd_delta_v[2] = (dd_delta * 0.0 + a * 0.0) + delta * 0.0;
  tau[2] = d1 / norm_d_delta_v;
  for (r1 = 0; r1 < 3; r1++) {
    d_tau_tmp[3 * r1] = d_delta_v[0] * d_delta_v[r1];
    d_tau_tmp[3 * r1 + 1] = d_delta_v[1] * d_delta_v[r1];
    d_tau_tmp[3 * r1 + 2] = d1 * d_delta_v[r1];
  }

  d_delta = rt_powd_snf(norm_d_delta_v, 3.0);
  for (r1 = 0; r1 < 3; r1++) {
    b_dd_delta_v[3 * r1] = dd_delta_v[0] * d_delta_v[r1];
    b_dd_delta_v[3 * r1 + 1] = dd_delta_v[1] * d_delta_v[r1];
    d = (d_tau_tmp[r1] * dd_delta_v[0] + d_tau_tmp[r1 + 3] * dd_delta_v[1]) +
      d_tau_tmp[r1 + 6] * dd_delta_v[2];
    b_dd_delta_v[3 * r1 + 2] = dd_delta_v[2] * d_delta_v[r1];
    e_phi[r1] = d;
    d_tau[r1] = dd_delta_v[r1] / norm_d_delta_v - d / d_delta;
  }

  phi = 0.0;
  for (r1 = 0; r1 < 3; r1++) {
    d_e_phi[r1] = (b_dd_delta_v[r1] * dd_delta_v[0] + b_dd_delta_v[r1 + 3] *
                   dd_delta_v[1]) + b_dd_delta_v[r1 + 6] * dd_delta_v[2];
    phi += dd_delta_v[r1] * dd_delta_v[r1];
  }

  d_gamma = rt_powd_snf(norm_d_delta_v, 5.0);
  for (r1 = 0; r1 < 3; r1++) {
    b_dd_delta_v[3 * r1] = e_phi[0] * 3.0 * d_delta_v[r1];
    b_dd_delta_v[3 * r1 + 1] = e_phi[1] * 3.0 * d_delta_v[r1];
    b_dd_delta_v[3 * r1 + 2] = e_phi[2] * 3.0 * d_delta_v[r1];
    dd_tau[r1] = (ddd_delta_v[r1] / norm_d_delta_v - d_e_phi[r1] / d_delta) -
      ((d_e_phi[r1] + d_delta_v[r1] * phi) + ((d_tau_tmp[r1] * ddd_delta_v[0] +
         d_tau_tmp[r1 + 3] * ddd_delta_v[1]) + d_tau_tmp[r1 + 6] * ddd_delta_v[2]))
      / d_delta;
  }

  for (r1 = 0; r1 < 3; r1++) {
    dd_tau[r1] += ((b_dd_delta_v[r1] * dd_delta_v[0] + b_dd_delta_v[r1 + 3] *
                    dd_delta_v[1]) + b_dd_delta_v[r1 + 6] * dd_delta_v[2]) /
      d_gamma;
  }

  //     %% Normal vector:
  //     %% Rho :
  //     %% Need to find varphi = g(phi)
  phi = 0.0;
  d_gamma = 0.0;
  m_12 = 0.0;
  b_t = 0.0;
  for (r1 = 0; r1 < 3; r1++) {
    signed char i;
    i = iv[r1 + 3];
    d_e_phi[r1] = delta_v[r1] + 0.0094868329805051412 * ((static_cast<double>
      (iv[r1]) * tau[0] + static_cast<double>(i) * tau[1]) + 0.0 * tau[2]);
    e_phi[r1] = d_delta_v[r1] + 0.0094868329805051412 * ((static_cast<double>
      (iv[r1]) * d_tau[0] + static_cast<double>(i) * d_tau[1]) + 0.0 * d_tau[2]);
    ddd_delta_v[r1] = dd_delta_v[r1] + 0.0094868329805051412 * ((static_cast<
      double>(iv[r1]) * dd_tau[0] + static_cast<double>(i) * dd_tau[1]) + 0.0 *
      dd_tau[2]);
    d = 0.0094868329805051412 * tau[r1];
    phi += delta_v[r1] * static_cast<double>(b[r1]);
    d_gamma += d * static_cast<double>(b_b[r1]);
    m_12 += delta_v[r1] * static_cast<double>(b_b[r1]);
    b_t += d * static_cast<double>(b[r1]);
  }

  b_gamma = phi - d_gamma;
  d_delta = m_12 + b_t;

  // func_g = atan2(gamma, psi);
  dd_s = b_gamma * b_gamma + d_delta * d_delta;

  //     %% s
  scale = 3.3121686421112381E-170;

  //     %% s_f:
  d = 0.0094868329805051412 * d_tau[0];
  d_gamma = d * 0.0;
  m_12 = d;
  d = 0.0094868329805051412 * dd_tau[0];
  d_sf = d * 0.0;
  delta = d;
  phi = std::abs(e_phi[0]);
  if (phi > 3.3121686421112381E-170) {
    d_s_tmp = 1.0;
    scale = phi;
  } else {
    b_t = phi / 3.3121686421112381E-170;
    d_s_tmp = b_t * b_t;
  }

  d = 0.0094868329805051412 * d_tau[1];
  d_gamma += d;
  m_12 += d * 0.0;
  d = 0.0094868329805051412 * dd_tau[1];
  d_sf += d;
  delta += d * 0.0;
  phi = std::abs(e_phi[1]);
  if (phi > scale) {
    b_t = scale / phi;
    d_s_tmp = d_s_tmp * b_t * b_t + 1.0;
    scale = phi;
  } else {
    b_t = phi / scale;
    d_s_tmp += b_t * b_t;
  }

  d = 0.0094868329805051412 * d_tau[2];
  d_gamma += d * 0.0;
  m_12 += d * 0.0;
  d = 0.0094868329805051412 * dd_tau[2];
  d_sf += d * 0.0;
  delta += d * 0.0;
  phi = std::abs(e_phi[2]);
  if (phi > scale) {
    b_t = scale / phi;
    d_s_tmp = d_s_tmp * b_t * b_t + 1.0;
    scale = phi;
  } else {
    b_t = phi / scale;
    d_s_tmp += b_t * b_t;
  }

  d_gamma = ((d_delta_v[0] + d_delta_v[1] * 0.0) + d1 * 0.0) - d_gamma;
  phi = ((d_delta_v[0] * 0.0 + d_delta_v[1]) + d1 * 0.0) + m_12;
  dd_delta = d_gamma * d_delta - phi * b_gamma;
  b_t = dd_delta / dd_s;
  d_gamma = ((((dd_delta_v[0] + dd_delta_v[1] * 0.0) + dd_delta_v[2] * 0.0) -
              d_sf) * d_delta - (((dd_delta_v[0] * 0.0 + dd_delta_v[1]) +
    dd_delta_v[2] * 0.0) + delta) * b_gamma) / dd_s - dd_delta * (2.0 * b_gamma *
    d_gamma + 2.0 * d_delta * phi) / (dd_s * dd_s);
  d_s_tmp = scale * std::sqrt(d_s_tmp);
  delta = d_s_tmp / b_t;
  phi = rt_powd_snf(b_t, 3.0);
  dd_delta = b_t * b_t;
  dd_s = ((e_phi[0] * ddd_delta_v[0] + e_phi[1] * ddd_delta_v[1]) + e_phi[2] *
          ddd_delta_v[2]) / d_s_tmp / dd_delta - d_s_tmp * d_gamma / phi;
  d_sf = norm_d_delta_v / b_t;
  d_gamma = ((d_delta_v[0] * dd_delta_v[0] + d_delta_v[1] * dd_delta_v[1]) + d1 *
             dd_delta_v[2]) / norm_d_delta_v / b_t - norm_d_delta_v * d_gamma /
    phi;

  //     %% Kappa:
  d_delta = delta * delta;

  //     %% Rotations:
  phi = std::sin(q_hat[0]);
  b_t = std::cos(q_hat[0]);
  e_phi[2] = d_e_phi[0] * tau[1] - d_e_phi[1] * tau[0];
  m_12 = 0.003 * (delta * e_phi[2] - d_sf * 4.44556451612904E-7 /
                  0.0094868329805051412 / 0.003);
  M[2] = m_12;
  M[1] = m_12;
  M[3] = 0.003 * (d_delta + d_sf * d_sf * 4.44556451612904E-7 /
                  9.000000000000006E-5 / 0.003);
  b_gamma = 99.999999999999986 * (theta - q_hat[0]);
  d_s_tmp = 99.999999999999986 * (varphi - q_hat[1]);
  ddd_delta_v[0] = ddd_delta_v[0] / dd_delta / d_delta;
  ddd_delta_v[1] = ddd_delta_v[1] / dd_delta / d_delta;
  M[0] = 0.003 * ((((d_e_phi[0] * d_e_phi[0] + d_e_phi[1] * d_e_phi[1]) +
                    d_e_phi[2] * d_e_phi[2]) + 0.00014818548387096798) +
                  0.52705666666666662);
  b_dd_delta_v[0] = -phi;
  b_dd_delta_v[3] = -b_t;
  b_dd_delta_v[6] = 0.0;
  b_dd_delta_v[1] = b_t;
  b_dd_delta_v[4] = -std::sin(q_hat[0]);
  b_dd_delta_v[7] = 0.0;
  scale = 0.0;
  for (r1 = 0; r1 < 3; r1++) {
    scale += (0.0 * b_dd_delta_v[3 * r1] + 9.81 * b_dd_delta_v[3 * r1 + 1]) *
      d_e_phi[r1];
    b_dd_delta_v[3 * r1 + 2] = b_iv[r1];
  }

  b_dd_delta_v[0] = b_t;
  b_dd_delta_v[3] = -phi;
  b_dd_delta_v[6] = 0.0;
  b_dd_delta_v[1] = phi;
  b_dd_delta_v[4] = b_t;
  b_dd_delta_v[7] = 0.0;
  a = 0.0;
  for (r1 = 0; r1 < 3; r1++) {
    a += ((0.0 * b_dd_delta_v[3 * r1] + 9.81 * b_dd_delta_v[3 * r1 + 1]) + 0.0 *
          b_dd_delta_v[3 * r1 + 2]) * tau[r1];
  }

  B[0] = ((-(((0.003 * d_e_phi[0] * tau[0] + 0.003 * d_e_phi[1] * tau[1]) +
              0.003 * d_e_phi[2] * tau[2]) * delta * dq_hat[1]) * dq_hat[0] +
           -(0.003 * (((d_e_phi[0] * tau[0] + d_e_phi[1] * tau[1]) + d_e_phi[2] *
                       tau[2]) * delta * dq_hat[0] + ((dd_s * e_phi[2] - d_gamma
    * 4.44556451612904E-7 / 0.0094868329805051412 / 0.003) + d_delta * (d_e_phi
    [0] * ddd_delta_v[1] - d_e_phi[1] * ddd_delta_v[0])) * dq_hat[1])) * dq_hat
           [1]) - 0.003 * scale) + u;
  B[1] = (-(((-0.003 * d_e_phi[0] * tau[0] + -0.003 * d_e_phi[1] * tau[1]) +
             -0.003 * d_e_phi[2] * tau[2]) * delta * dq_hat[0]) * dq_hat[0] +
          -(0.003 * (delta * dd_s + d_sf * d_gamma * 4.44556451612904E-7 /
                     9.000000000000006E-5 / 0.003) * dq_hat[1]) * dq_hat[1]) -
    0.003 * (a * delta);
  if (std::abs(m_12) > std::abs(M[0])) {
    r1 = 1;
    r2 = 0;
  } else {
    r1 = 0;
    r2 = 1;
  }

  phi = M[r2] / M[r1];
  dd_delta = M[r1 + 2];
  m_12 = (B[r2] - B[r1] * phi) / (M[r2 + 2] - phi * dd_delta);
  b_t = t - t_prev;

  //  + ;
  d = q_hat[0];
  d1 = q_hat[0];
  phi = q_hat[1];
  d_gamma = q_hat[1];
  q_hat[0] = d + (dq_hat[0] + 10.0 * (theta - d1)) * b_t;
  q_hat[1] = phi + (dq_hat[1] + 10.0 * (varphi - d_gamma)) * b_t;
  t_prev = t;
  d = dq_hat[0] + ((B[r1] - m_12 * dd_delta) / M[r1] + b_gamma) * b_t;
  dq_hat[0] = d;
  dq[0] = d;
  q[0] = q_hat[0];
  d = dq_hat[1] + (m_12 + d_s_tmp) * b_t;
  dq_hat[1] = d;
  dq[1] = d;
  q[1] = q_hat[1];
}

void observer_init()
{
  t_prev = 0.0;
  dq_hat[0] = 0.0;
  dq_hat[1] = 0.0;
}

void spline_phi_not_empty_init()
{
  spline_phi_not_empty = false;
}

// End of code generation (observer.cpp)
