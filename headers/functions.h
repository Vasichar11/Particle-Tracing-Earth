#pragma once
#include <cmath>
#include <iostream>
#include <tuple>
#include "constants.h"
#include "common.h"


//-----------For Bell-----------//

//Stix tuple, needs wps,wc and returns SDPRL.
std::tuple<real, real, real, real, real> stix_parameters(real wce, real wcO, real wcH, real wcHe, real wpse, real wpsO, real wpsH, real wpsHe);
    //returns tuple(S, D, P, R, L).

//Magnetic dipole field
real Bmag_dipole(real latitude);
    //returns Bmag.

//Solve dispersion relation.
std::tuple<real, real, real, real> dispersion(real S, real P, real R, real L, real D);
    //returns tuple(mu, kappa, kx, kz). 

//Estimate dwh_ds
real dwh_dsf(real w_h, real latitude);
    //returns dwh_ds.

//Estimate Bell parameters.
void Bell_params(const real ppar,const real pper, const real Bxw, const real Byw, const real Exw, const real Eyw, const real Ezw, const real kz, const real kx, const real w_h, real gama, real &w1, real &w2, real &wtau_sq, real &R1, real &R2, real &beta);
    //returns by reference  gama,  w1,  w2,  wtau_sq,  R1,  R2,  beta.

//Estimate resonant velocity
void vres_f(const real kz, const real w_h, const real alpha, real &v_para_res, real &E_res);
    //returns by reference v_para_res,  E_res.

//Compute field components.
void whistlers(int64_t p, int64_t i, real mu, real P, real D, real S, real kz, real &Bxw,real &Byw,real &Bzw, real &Exw,real &Eyw,real &Ezw);
    //returns by reference  Bxw, Byw, Bzw,  Exw, Eyw, Ezw.



//-------For Ray Tracing-------//

//Return Factors, aeq of runge kutta and kz vector, only when particle in packet.
void f_packet (real &Fpar, real &Fper, real &Ftheta, real &aeq_rk, real &kz, const real pper_tmp, const real ppar_tmp, const real eta_tmp, const real aeqsu_tmp, const real alpha_tmp, const real gama, const real w_h, const real p_mag, const real kx_ray, const real kz_ray, const real kappa_ray, const real Bw_ray, const real Bzw, const real Ezw, const real w1, const real w2, const real R1, const real R2 );


//Returns quantities needed, regardless WPI.
void f_always(real &p_mag, real &gama, real &w_h, real &dwh_ds, const real latitude_tmp, const real ppar_tmp, const real pper_tmp);


//Estimate resonant velocity 
real vres_f(real w_wave, real kz, real w_h, real alpha); 