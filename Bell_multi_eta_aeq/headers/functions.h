#pragma once
#include <cmath>
#include <iostream>
#include <tuple>
#include "constants.h"
#include "common.h"


//Prototypes of functions.

//Stix tuple, needs wps,wc and returns SDPRL.
std::tuple<real, real, real, real, real> stix_parameters(real w, real wce, real wcO, real wcH, real wcHe, real wpse, real wpsO, real wpsH, real wpsHe);
    //returns tuple(S, D, P, R, L).

//Magnetic dipole field
real Bmag_dipole(real lamda);
    //returns Bmag.

//Solve dispersion relation.
std::tuple<real, real, real, real> dispersion(real S, real P, real R, real L, real D, real th, real w_wave);
    //returns tuple(mu, kappa, kx, kz). 

//Estimate dwh_ds
real dwh_dsf(real w_h, real lamda);
    //returns dwh_ds.

//Estimate Bell parameters.
void Bell_params(const real m_e, const real q_e, const real ppar, const real pper, const real Bxw, const real Byw, const real Exw, const real Eyw, const real Ezw, const real kz, const real kx, const real w_h, real &gama, real &w1, real &w2, real &wtau_sq, real &R1, real &R2, real &beta);
    //returns by reference  gama,  w1,  w2,  wtau_sq,  R1,  R2,  beta.

//Estimate resonant velocity
void vres_f(const real w, const real kz, const real w_h, const real alpha, real &v_para_res, real &E_res);
    //returns by reference v_para_res,  E_res.

//Compute field components.
void whistlers(int64_t p, int64_t i, real mu, real P, real D, real S, real kz, real zeta, real time, real &Bxw,real &Byw,real &Bzw, real &Exw,real &Eyw,real &Ezw);
    //returns by reference  Bxw, Byw, Bzw,  Exw, Eyw, Ezw.


