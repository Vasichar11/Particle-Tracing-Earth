#pragma once
#include "common.h"
#include "constants.h"
#include <cmath>

// Function to calculate all needed slopes for the particular step (adiabatic motion).
void slopes(real &k, real &l, real &m, real &o, real &p, const real ppar_tmp, const real pper_tmp, const real latitude_tmp, const real w_h, const real dwh_ds, const real gama);

// Overloaded Function to calculate all needed slopes for the particular step (overloaded for WPI).
void slopes(real &k, real &l, real &m, real &n, real &o, real &p, real &q, const real ppar_tmp, const real pper_tmp, const real latitude_tmp, const real eta_tmp, const real alpha_tmp, const real aeq_tmp, const real p_mag, const real w_h, const real dwh_ds, const real gama, const real kz, const real kappa, const real wtau_sq, const real w1, const real w2, const real R1, const real R2, const real beta, const real Bw_out);

// Overloaded Function to calculate all needed slopes for the particular step (overloaded for WPI in Ray Tracing).
void slopes(real &k, real &l, real &m, real &n, real &o, real &p, real &q,  const real ppar_tmp, const real pper_tmp, const real alpha_tmp, const real latitude_tmp, const real eta_tmp, const real Fpar, const real Fper, const real Ftheta, const real gama, const real w_h, const real dwh_ds, const real kz, real Bw_ray, real p_mag, real aeqsu_tmp, real kappa_ray);
