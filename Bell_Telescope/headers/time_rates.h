#pragma once
#include "common.h"
#include "constants.h"
#include <cmath>

//Function to calculate all needed slopes for the particular step(adiabatic motion).
void slopes(real &k, real &l, real &m, real &o, real &p, real ppar_tmp, real pper_tmp, real lamda_tmp, real w_h, real dwh_ds, real gama);

//Function to calculate all needed slopes for the particular step(overloaded for WPI).
void slopes(real &k, real &l, real &m, real &n, real &o, real &p, real &q, real ppar_tmp, real pper_tmp, real lamda_tmp, real eta_tmp, real alpha_tmp, real aeq_tmp, real p_mag, real w_h, real dwh_ds, real gama, real kz, real kappa, real wtau_sq, real w1, real w2, real R1, real R2, real beta, real Bw_out);
