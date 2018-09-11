#ifndef __COOK_TORRANCE_H__
#define __COOK_TORRANCE_H__

#include "fresnel.h"
#include "cook_torrance_common.h"


inline double BeckmannDistribution(K_DBL mu, K_DBL phi, const bool phi_is_cos_phi = false){
  K_DBL cos_phi = phi_is_cos_phi ? phi : std::cos(phi),
        mu_sqrd = mu*mu;
  K_DBL cos_phi_sqrd = cos_phi*cos_phi;

  return ( std::exp( -(1.0-cos_phi_sqrd)/(cos_phi_sqrd*mu_sqrd) ) / (M_PI*mu_sqrd*cos_phi_sqrd*cos_phi_sqrd) );
}

inline double BeckmannDistribution(K_DBL &mu, const vec3d_t &H, const vec3d_t &N){
  K_DBL cos_phi = DOT_MACRO(H, N);
  return BeckmannDistribution(mu, cos_phi, true);
}

//inline double BeckmannDistribution_d_dmu(K_DBL &mu, K_DBL &phi, const bool phi_is_cos_phi = false){
//  K_DBL cos_phi = phi_is_cos_phi ? phi : std::cos(phi),
//               tan_phi = phi_is_cos_phi ? std::tan( std::acos(phi) ) : std::tan(phi),
//               mu_sqrd = mu*mu;
//  K_DBL tan_phi_sqrd = tan_phi*tan_phi;
//  return ( -2.0*(mu_sqrd-tan_phi_sqrd)*std::exp(-tan_phi_sqrd / mu_sqrd) ) /
//         (std::pow(cos_phi, 4)*mu_sqrd*mu_sqrd*mu);
//}

inline double BeckmannDistribution_d_dmu(K_DBL &mu, K_DBL &phi, const bool phi_is_cos_phi = false){
  K_DBL cos_phi = phi_is_cos_phi ? phi : std::cos(phi),
        mu_sqrd = mu*mu;
  K_DBL cos_phi_sqrd = cos_phi*cos_phi;
  return ( -2.0 * (cos_phi_sqrd*(mu_sqrd+1.0) - 1.0) * std::exp( -(1.0-cos_phi_sqrd)/(cos_phi_sqrd*mu_sqrd) ) ) /
         (M_PI * cos_phi_sqrd*cos_phi_sqrd*cos_phi_sqrd * mu_sqrd*mu_sqrd*mu);
}

inline double BeckmannDistribution_d_dmu(K_DBL &mu, const vec3d_t &H, const vec3d_t &N){
  K_DBL cos_phi = DOT_MACRO(H, N);
  return BeckmannDistribution_d_dmu(mu, cos_phi, true);
}


inline double GeometricAttenuation(const vec3d_t &H, const vec3d_t &L, const vec3d_t &N, const vec3d_t &V){
  K_DBL temp = ( 2.0 * DOT_MACRO(H, N) ) / DOT_MACRO(H, V);
  return std::min( 1.0, std::min( temp*DOT_MACRO(N, V), temp*DOT_MACRO(N, L) ) );
}


inline double CookTorrance(const vec3d_t &H, const vec3d_t &L, const vec3d_t &N, const vec3d_t &V,
                           K_DBL &n, K_DBL &k, K_DBL &mu){
  K_DBL D = BeckmannDistribution(mu, H, N),
        G = GeometricAttenuation(N, H, V, L),
        F = FresnelSimple(n, k, H, V);
  return (D*F*G) / ( 4.0*DOT_MACRO(N, L)*DOT_MACRO(N, V) );
}


inline void CookTorrance_deriv(const vec3d_t &H, const vec3d_t &L, const vec3d_t &N, const vec3d_t &V,
                               K_DBL &n, K_DBL &k, K_DBL &mu,
                               double &dR_dn, double &dR_dk, double &dR_dmu){
  K_DBL D = BeckmannDistribution(mu, H, N),
        G = GeometricAttenuation(N, H, V, L),
        F = FresnelSimple(n, k, H, V);
  K_DBL bot = 1.0 / ( 4.0*DOT_MACRO(N, L)*DOT_MACRO(N, V) );

  double dF_dn, dF_dk;
  FresnelSimple_deriv(n, k, H, V, dF_dn, dF_dk);
  K_DBL temp = (D*G) * bot;
  dR_dn = temp * dF_dn;
  dR_dk = temp * dF_dk;

  K_DBL dD_dmu = BeckmannDistribution_d_dmu(mu, H, N);
  dR_dmu = (F*G*dD_dmu) * bot;
}

#endif //__COOK_TORRANCE_H__

