#pragma once

#include <cmath>
#include <iostream>
#include "common.h"
#include "constants.h"

//Return Factors, aeq of runge kutta and kz vector, only when particle in packet.
void f_packet (real &Fpar, real &Fper, real &Ftheta, real &aeq_rk, real &kz, const real pper_tmp, const real ppar_tmp, const real eta_tmp, const real aeqsu_tmp, const real alpha_tmp, const real gama, const real w_h, const real p_mag, const real kx_ray, const real kz_ray, const real kappa_ray, const real Bw_ray, const real Bzw, const real Ezw, const real w1, const real w2, const real R1, const real R2 );


//Returns quantities needed, regardless WPI.
void f_always(real &p_mag, real &gama, real &w_h, real &dwh_ds, const real lamda_tmp, const real ppar_tmp, const real pper_tmp);


//Magnetic dipole field 
real Bmag_dipole(real lamda);                                     


//Estimate resonant velocity 
real vres_f(real w_wave, real kz, real w_h, real alpha); 
