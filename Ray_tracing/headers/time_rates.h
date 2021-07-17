#pragma once

#include "common.h"
#include "constants.h"
#include <cmath>


//Returns slopes of Runge kutta in each step.
void slopes(real &k, real &l, real &m, real &n, real &o, real &p, const real ppar_tmp, const real pper_tmp, const real alpha_tmp, const real lamda_tmp, const real eta_tmp, const real Fpar, const real Fper, const real Ftheta, const real gama, const real w_h, const real dwh_ds, const real kz);
