#ifndef __FRESNEL_H__
#define __FRESNEL_H__

#include <cmath>
#include "cook_torrance_common.h"


// Used in Ianiro paper
// From book: Siegel and Howell - 3.2.2.1 Electromagnetic Relations for Incidence on an Absorbing Medium
// All equation references are from Siegel and Howell, 5th Edition
inline double FresnelComplex(K_DBL &n, K_DBL &k, K_DBL &theta){
  K_DBL cos_t = std::cos(theta),
        sin_t = std::sin(theta);
  K_DBL n_sq = n*n,
        k_sq = k*k;
  K_DBL n_sq_add_k_sq = n_sq+k_sq;
  K_DBL n_sq_add_k_sq_inv = 1.0/n_sq_add_k_sq;

  //TODO - when n or k get too low, the operands of the square roots becomes negative and we get -nan outputs.
  K_DBL tmp0 = sin_t*sin_t * n_sq_add_k_sq_inv;
  K_DBL alpha = std::sqrt( (1.0 + tmp0) - (4.0*n_sq * n_sq_add_k_sq_inv) * tmp0 ); //ref (3.14c)
  K_DBL tmp1 = (n_sq-k_sq) * n_sq_add_k_sq_inv,
        tmp2 = n_sq_add_k_sq/n_sq;
  K_DBL beta = std::sqrt( (tmp2*0.5) * (tmp1 - tmp0 + alpha) );  //ref (3.14d)
  K_DBL gamma = (tmp1*beta) + (2.0*n*k * n_sq_add_k_sq_inv) * std::sqrt(tmp2*alpha - beta);  //ref (3.14e)

  K_DBL n_beta = n*beta;
  K_DBL n_beta_min_cos_t = n_beta - cos_t,
        tmp3 = n_sq_add_k_sq*alpha;
  K_DBL tmp4 = tmp3 - n_sq*beta*beta,
        n_beta_add_cos_t = n_beta + cos_t;
  K_DBL rho_perp = (n_beta_min_cos_t*n_beta_min_cos_t + tmp4) /
                   (n_beta_add_cos_t*n_beta_add_cos_t + tmp4); //ref (3.14a)

  K_DBL n_gamma = n*gamma,
        alpha_div_cos_t = alpha/cos_t;
  K_DBL tmp5 = n_gamma - alpha_div_cos_t,
        tmp6 = tmp3 - n_sq*gamma*gamma,
        tmp7 = n_gamma + alpha_div_cos_t;
  K_DBL rho_para = (tmp5*tmp5 + tmp6) /
                   (tmp7*tmp7 + tmp6);  //ref (3.14b)

  return ( (rho_perp+rho_para) * 0.5 );
}

inline double FresnelEmissivityComplex(K_DBL &n, K_DBL &k, K_DBL &theta){
  return ( 1.0 - FresnelComplex(n, k, theta) );
}


// From book: Siegel and Howell - 3.2.2.2 Reflectivity and Emissivity Relations for Metals (Large k)
// All equation references are from Siegel and Howell, 5th Edition
inline double FresnelSimple(K_DBL &n, K_DBL &k, K_DBL &theta, const bool theta_is_cos_theta = false){
  K_DBL z = theta_is_cos_theta ? theta : std::cos(theta); //z = cos(theta)
  K_DBL nz_min_one = n*z - 1.0,
        nz_add_one = n*z + 1.0,
        k_sqrd = k*k,
        n_min_z = n-z,
        n_add_z = n+z;
  K_DBL kz_sqrd = k_sqrd*z*z;

  K_DBL rho_para = (nz_min_one*nz_min_one + kz_sqrd) /
                   (nz_add_one*nz_add_one + kz_sqrd); //ref (3.15a)
  K_DBL rho_perp = (n_min_z*n_min_z + k_sqrd) /
                   (n_add_z*n_add_z + k_sqrd); //ref (3.15b)

  return ( (rho_para+rho_perp) * 0.5 ); //ref (3.16)
}

inline double FresnelSimple(K_DBL &n, K_DBL &k, const vec3d_t &H, const vec3d_t &V){
  K_DBL cos_theta = DOT_MACRO(H, V);
  return FresnelSimple(n, k, cos_theta, true);
}

inline double FresnelEmissivitySimple(K_DBL &n, K_DBL &k, K_DBL &theta){
  return ( 1.0 - FresnelSimple(n, k, theta) );
}


inline void FresnelSimple_deriv(K_DBL &n, K_DBL &k, K_DBL &theta, double &d_dn, double &d_dk,
                                const bool theta_is_cos_theta = false){
  K_DBL z = theta_is_cos_theta ? theta : std::cos(theta); //z = cos(theta)
  K_DBL z_sqrd = z*z,
        k_sqrd = k*k,
        nz_min_one = n*z - 1.0,
        nz_add_one = n*z + 1.0,
        n_min_z = n-z,
        n_add_z = n+z;
  K_DBL kz_sqrd = k_sqrd*z_sqrd;

  double val[4];
  K_DBL com0 = 1.0 / (nz_add_one*nz_add_one + kz_sqrd);
  K_DBL com0_sqrd = com0*com0;
  val[0] = (z*nz_min_one) * com0;
  K_DBL com1 = nz_min_one*nz_min_one + kz_sqrd;
  val[1] = ( (z*nz_add_one) * com1 ) * com0_sqrd;
  K_DBL com2 = 1.0 / (n_add_z*n_add_z + k_sqrd);
  K_DBL com2_sqrd = com2*com2;
  val[2] = n_min_z * com2;
  K_DBL com3 = n_min_z*n_min_z + k_sqrd;
  val[3] = (n_add_z * com3) * com2_sqrd;
  d_dn = val[0] - val[1] + val[2] - val[3];
  
  K_DBL k_z_sqrd = z_sqrd*k;
  val[0] = k_z_sqrd * com0;
  val[1] = (k_z_sqrd * com1) * com0_sqrd;
  val[2] = k * com2;
  val[3] = (k * com3) * com2_sqrd;
  d_dk = val[0] - val[1] + val[2] - val[3];
}

inline void FresnelSimple_deriv(K_DBL &n, K_DBL &k, const vec3d_t &H, const vec3d_t &V, double &d_dn, double &d_dk,
                                const bool theta_is_cos_theta = false){
  K_DBL cos_theta = DOT_MACRO(H, V);
  FresnelSimple_deriv(n, k, cos_theta, d_dn, d_dk, true);
}


// From book: Siegel and Howell - 3.2.2.2 Reflectivity and Emissivity Relations for Metals (Large k)
// All equation references are from Siegel and Howell, 5th Edition
inline double FresnelEmissivitySimpleDirect(K_DBL &n, K_DBL &k, K_DBL &theta, const bool theta_is_cos_theta = false){
  K_DBL cos_t = theta_is_cos_theta ? theta : std::cos(theta);

  K_DBL emis_para = (4.0*n*cos_t) /
                    (cos_t*cos_t + 2.0*n*cos_t + n*n + k*k); //ref (3.17a)
  K_DBL emis_perp = (4.0*n*cos_t) /
                    ( (n*n + k*k) * cos_t*cos_t + 2.0*n*cos_t + 1.0 ); //ref (3.17b)

  return ( (emis_perp+emis_para) * 0.5 ); //ref (3.18)
}


// From paper: Fresnel Term Approximations for Metals - Istvan Lazanyi, Laszlo Szimay-Kalos
inline double FresnelMetalApprox(K_DBL &n, K_DBL &k, K_DBL &theta, K_DBL &err_a, K_DBL &err_gamma,
                                 const bool theta_is_cos_theta = false){
  K_DBL cos_t = theta_is_cos_theta ? theta : std::cos(theta),
                k_sqrd = k*k;
  K_DBL top = (n-1.0)*(n-1.0) + 4.0*n*std::pow(1.0-cos_t, 5) + k_sqrd,
               bot = (n+1.0)*(n+1.0) + k_sqrd,
               err = err_a*cos_t*std::pow(1.0-cos_t, err_gamma);
  return top/bot - err;
}

inline double FresnelMetalApprox(K_DBL &n, K_DBL &k, const vec3d_t &H, const vec3d_t &V, K_DBL &err_a, K_DBL &err_gamma){
  K_DBL cos_theta = DOT_MACRO(H, V);
  return FresnelMetalApprox(n, k, cos_theta, err_a, err_gamma, true);
}

inline void FresnelMetalApprox_deriv(K_DBL &n, K_DBL &k, K_DBL &theta, double &d_dn, double &d_dk,
                                     const bool theta_is_cos_theta = false){
  K_DBL cos_t = theta_is_cos_theta ? theta : std::cos(theta);
  K_DBL cos_t_sqrd = cos_t*cos_t,
                     k_sqrd = k*k,
                     n_plus = n+1.0;

  K_DBL bot = std::pow(n_plus*n_plus + k_sqrd, 2);
  K_DBL common = ( 4.0*cos_t*(cos_t_sqrd*cos_t_sqrd - 5.0*cos_t_sqrd*cos_t + 10.0*cos_t*(cos_t-1.0) + 5.0) ) / bot;
               
  d_dn = common*(n*n - k_sqrd - 1.0);
  d_dk = 2.0*n*common*k;
}

inline void FresnelMetalApprox_deriv(K_DBL &n, K_DBL &k, const vec3d_t &H, const vec3d_t &V, double &d_dn, double &d_dk,
                                     const bool theta_is_cos_theta = false){
  K_DBL cos_theta = DOT_MACRO(H, V);
  FresnelMetalApprox_deriv(n, k, cos_theta, d_dn, d_dk, true);
}

#endif //__FRESNEL_H__

