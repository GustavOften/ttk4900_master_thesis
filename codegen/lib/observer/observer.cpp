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

//
//  function [q,dq] = observer(t,theta,varphi,u)

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

  // 'observer:2' epsilon = 0.01;
  // 'observer:3' k1 = 2;
  // 'observer:3' k2 = 1;
  // 'observer:5' if isempty(spline_phi)
  if (!spline_phi_not_empty) {
    // 'observer:6' t_prev = 0;
    // 'observer:7' q_hat = [theta;varphi];
    q_hat[0] = theta;
    q_hat[1] = varphi;

    // 'observer:8' dq_hat = [0;2];
    // 'observer:9' solutions = coder.load('solutions_riccati_and_dphi');
    // 'observer:10' [breaks,coefs,~,~,dim] = unmkpp(solutions.varphi_to_phi);
    // 'observer:10' ~
    // 'observer:10' ~
    // 'observer:11' spline_phi = mkpp(breaks,coefs,dim);
    std::memcpy(&spline_phi.breaks[0], &dv[0], 500U * sizeof(double));
    std::memcpy(&spline_phi.coefs[0], &dv1[0], 1996U * sizeof(double));
    spline_phi_not_empty = true;
  }

  // n_varphi = solutions.n_varphi;
  // 'observer:14' q = q_hat;
  // 'observer:15' dq = dq_hat;
  // 'observer:17' a = 0.1095;
  // Constant for describing the butterfly frame
  // 'observer:18' b = 0.0405;
  // Constant for describing the butterfly frame
  // 'observer:19' g = 9.81;
  // Gravity
  // 'observer:20' J_s = 4.44556451612904e-07;
  // 5.48e-7;%5.48e-7; %Mass moment of inertia of the ball
  // 'observer:21' J_f = 1.58117e-3;
  // Mass moment of inertia of the frame, from article.
  // 'observer:22' m_b = 3.0e-3;
  // Mass of the ball
  // R_b = sqrt(16.55e-3^2-12.5e-3^2); From paper.
  // 'observer:24' R_b = sqrt(16.5e-3^2-13.5e-3^2);
  // Measured on the robot
  // 'observer:26' phi = ppval(spline_phi,mod(q(2),2*pi));
  phi = ppval(spline_phi.breaks, spline_phi.coefs, b_mod(q_hat[1]));

  // 'observer:27' phi = phi(1,1);
  // Need to do this becaouse the simulink coder is strange...
  //     %% e_phi
  // 'observer:30' e_phi = [sin(phi);cos(phi);0];
  dd_s = std::sin(phi);
  d_sf = std::cos(phi);

  // 'observer:31' d_e_phi = [cos(phi);-sin(phi);0];
  // 'observer:32' dd_e_phi = [-sin(phi);-cos(phi);0];
  // 'observer:33' ddd_e_phi = [-cos(phi);sin(phi);0];
  //     %% Delta
  // 'observer:36' delta = (a-b*cos(2*phi));
  d_gamma = std::cos(2.0 * phi);
  delta = 0.1095 - 0.0405 * d_gamma;

  // 'observer:37' delta_v = delta*e_phi;
  // 'observer:38' d_delta = b*2*sin(2*phi);
  phi = std::sin(2.0 * phi);
  d_delta = 0.081 * phi;

  // 'observer:39' dd_delta = b*4*cos(2*phi);
  dd_delta = 0.162 * d_gamma;

  // 'observer:40' ddd_delta = -b*8*sin(2*phi);
  scale = -0.324 * phi;

  // 'observer:41' d_delta_v = d_delta*e_phi+delta*d_e_phi;
  // 'observer:42' dd_delta_v = dd_delta*e_phi+2*d_delta*d_e_phi+delta*dd_e_phi; 
  a = 2.0 * d_delta;

  // 'observer:43' ddd_delta_v = ddd_delta*e_phi+ 3*dd_delta*d_e_phi + ...
  // 'observer:44'         3*d_delta*dd_e_phi + delta*ddd_e_phi;
  phi = 3.0 * dd_delta;
  d_gamma = 3.0 * d_delta;
  b_t = delta * -d_sf;
  ddd_delta_v[0] = ((scale * dd_s + phi * d_sf) + d_gamma * -dd_s) + b_t;
  m_12 = delta * dd_s;
  ddd_delta_v[1] = ((scale * d_sf + phi * -dd_s) + d_gamma * -d_sf) + m_12;
  ddd_delta_v[2] = ((scale * 0.0 + phi * 0.0) + d_gamma * 0.0) + delta * 0.0;

  // 'observer:45' norm_d_delta_v = sqrt(d_delta^2+delta^2);
  norm_d_delta_v = std::sqrt(d_delta * d_delta + delta * delta);

  //     %% Tau
  // 'observer:48' tau = d_delta_v/norm_d_delta_v;
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

  // 'observer:49' d_tau = dd_delta_v/norm_d_delta_v - d_delta_v*d_delta_v'*dd_delta_v/norm_d_delta_v^3; 
  for (r1 = 0; r1 < 3; r1++) {
    d_tau_tmp[3 * r1] = d_delta_v[0] * d_delta_v[r1];
    d_tau_tmp[3 * r1 + 1] = d_delta_v[1] * d_delta_v[r1];
    d_tau_tmp[3 * r1 + 2] = d1 * d_delta_v[r1];
  }

  d_delta = rt_powd_snf(norm_d_delta_v, 3.0);

  // 'observer:50' dd_tau = ddd_delta_v/norm_d_delta_v - dd_delta_v*d_delta_v'*dd_delta_v/norm_d_delta_v^3 - ... 
  // 'observer:51'         (dd_delta_v*d_delta_v'*dd_delta_v+d_delta_v*(dd_delta_v'*dd_delta_v)+d_delta_v*d_delta_v'*ddd_delta_v)/... 
  // 'observer:52'         norm_d_delta_v^3 + d_delta_v*d_delta_v'*dd_delta_v*3*d_delta_v'*dd_delta_v/norm_d_delta_v^5; 
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
  // 'observer:54' normal = [0 -1 0;1 0 0;0 0 0]*tau;
  // 'observer:55' d_normal = [0 -1 0;1 0 0;0 0 0]*d_tau;
  // 'observer:56' dd_normal = [0 -1 0;1 0 0;0 0 0]*dd_tau;
  //     %% Rho :
  // 'observer:60' rho = delta_v + R_b*normal;
  // 'observer:61' d_rho = d_delta_v + R_b*d_normal;
  // 'observer:62' dd_rho = dd_delta_v + R_b*dd_normal;
  //     %% Need to find varphi = g(phi)
  // 'observer:65' gamma = delta_v'*[1;0;0]-R_b*tau'*[0;1;0];
  phi = 0.0;
  d_gamma = 0.0;

  // 'observer:66' psi = delta_v'*[0;1;0]+R_b*tau'*[1;0;0];
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

  // 'observer:67' d_gamma = d_delta_v'*[1;0;0]-R_b*d_tau'*[0;1;0];
  // 'observer:68' d_psi = d_delta_v'*[0;1;0]+R_b*d_tau'*[1;0;0];
  // 'observer:69' dd_gamma = dd_delta_v'*[1;0;0]-R_b*dd_tau'*[0;1;0];
  // 'observer:70' dd_psi = dd_delta_v'*[0;1;0]+R_b*dd_tau'*[1;0;0];
  // func_g = atan2(gamma, psi);
  // 'observer:72' d_func_g = (d_gamma*psi-d_psi*gamma)/(gamma^2+psi^2);
  // 'observer:73' dd_func_g = (dd_gamma*psi-dd_psi*gamma)/(psi^2+gamma^2)-(d_gamma*psi-d_psi*gamma)*(2*gamma*d_gamma+2*psi*d_psi)/(gamma^2+psi^2)^2; 
  dd_s = b_gamma * b_gamma + d_delta * d_delta;

  //     %% s
  // 'observer:76' d_s = norm(d_rho)/d_func_g;
  scale = 3.3121686421112381E-170;

  // 'observer:77' dd_s = d_rho'*dd_rho/norm(d_rho)/d_func_g^2 - norm(d_rho)*dd_func_g/d_func_g^3; 
  //     %% s_f:
  // 'observer:80' d_sf = norm_d_delta_v/d_func_g;
  // 'observer:81' dd_sf = d_delta_v'*dd_delta_v/norm_d_delta_v/d_func_g - norm_d_delta_v*dd_func_g/d_func_g^3; 
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
  // 'observer:84' kappa = dd_rho/d_func_g^2/d_s^2;
  d_delta = delta * delta;

  //     %% Rotations:
  // 'observer:87' R = [cos(q(1)) -sin(q(1)) 0;sin(q(1)) cos(q(1)) 0;0 0 1];
  phi = std::sin(q_hat[0]);
  b_t = std::cos(q_hat[0]);

  // 'observer:88' dR = [-sin(q(1)) -cos(q(1)) 0;cos(q(1)) -sin(q(1)) 0;0 0 0];
  // 'observer:90' rhoxtau = cross(rho,tau);
  e_phi[2] = d_e_phi[0] * tau[1] - d_e_phi[1] * tau[0];

  // 'observer:91' rhoxkappa = cross(rho,kappa);
  // 'observer:93' m_11 = m_b*(rho'*rho+J_s/m_b+J_f/m_b);
  // 'observer:94' m_12 = m_b*(d_s*rhoxtau(3)-d_sf*J_s/R_b/m_b);
  m_12 = 0.003 * (delta * e_phi[2] - d_sf * 4.44556451612904E-7 /
                  0.0094868329805051412 / 0.003);

  // 'observer:95' m_22 = m_b*(d_s^2+d_sf^2*J_s/R_b^2/m_b);
  // 'observer:96' M = [m_11 m_12;m_12 m_22];
  M[2] = m_12;
  M[1] = m_12;
  M[3] = 0.003 * (d_delta + d_sf * d_sf * 4.44556451612904E-7 /
                  9.000000000000006E-5 / 0.003);

  // 'observer:98' c_11 = m_b*rho'*tau*d_s*dq(2);
  // 'observer:99' c_12 = m_b*(rho'*tau*d_s*dq(1)+(dd_s*rhoxtau(3)-dd_sf*J_s/R_b/m_b+d_s^2*rhoxkappa(3))*dq(2)); 
  // 'observer:100' c_21 = -m_b*rho'*tau*d_s*dq(1);
  // 'observer:101' c_22 = m_b*(d_s*dd_s+d_sf*dd_sf*J_s/R_b^2/m_b)*dq(2);
  // 'observer:102' C = [c_11 c_12;c_21 c_22];
  // 'observer:105' G = m_b*[[0 g 0]*dR*rho;[0 g 0]*R*tau*d_s];
  // 'observer:106' ddq_hat_pri = M\(-C*dq-G+[u;0])+k1/epsilon^2*([theta;varphi]-q_hat); 
  b_gamma = 20000.0 * (theta - q_hat[0]);
  d_s_tmp = 20000.0 * (varphi - q_hat[1]);
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

  // 'observer:107' dq_hat_post = dq_hat + ddq_hat_pri*(t-t_prev);
  b_t = t - t_prev;

  //  + ;
  // 'observer:108' q_hat_post = q_hat + (dq_hat + k2/epsilon*([theta;varphi]-q_hat))*(t-t_prev); 
  d = q_hat[0];
  d1 = q_hat[0];
  phi = q_hat[1];
  d_gamma = q_hat[1];
  q_hat[0] = d + (dq_hat[0] + 100.0 * (theta - d1)) * b_t;
  q_hat[1] = phi + (dq_hat[1] + 100.0 * (varphi - d_gamma)) * b_t;

  // 'observer:109' t_prev = t;
  t_prev = t;

  // 'observer:110' dq_hat = dq_hat_post;
  // 'observer:111' q_hat = q_hat_post;
  // 'observer:112' dq = dq_hat;
  // 'observer:113' q = q_hat;
  d = dq_hat[0] + ((B[r1] - m_12 * dd_delta) / M[r1] + b_gamma) * b_t;
  dq_hat[0] = d;
  dq[0] = d;
  q[0] = q_hat[0];
  d = dq_hat[1] + (m_12 + d_s_tmp) * b_t;
  dq_hat[1] = d;
  dq[1] = d;
  q[1] = q_hat[1];
}

//
//  function [q,dq] = observer(t,theta,varphi,u)

void observer_init()
{
  t_prev = 0.0;
  dq_hat[0] = 0.0;
  dq_hat[1] = 2.0;
}

void spline_phi_not_empty_init()
{
  spline_phi_not_empty = false;
}

// End of code generation (observer.cpp)
