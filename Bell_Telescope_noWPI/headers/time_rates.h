#pragma once
#include "common.h"
#include "constants.h"
#include <cmath>


//Time rates functions prototypes.

//Equation 1: parallel(to B) speed uz=dz/dt, variation with time
real z_rk(real ppar_tmp, real m_e, real gama);			

//Equation 2: parallel momentum dpz/dt, variation with time
real p_par_rk(real pper_tmp,real eta_tmp, real m_e, real kz, real w_h, real dwh_ds, real gama, real wtau_sq);

//Equation 3: differentiate perpendicular momentum, variation with time
real p_per_rk(real ppar_tmp,real pper_tmp,real eta_tmp, real m_e, real w_h, real dwh_ds, real gama, real w1, real w2, real R1, real R2, real beta);

//Equation 4: eta: angle between BwR and u_per, variation with time
real eta_rk(real m_e,real w_wave,real kz,real ppar_tmp, real w_h, real gama);

//Equation 5: lamda variation with time
real lamda_rk(real ppar_tmp,real lamda_tmp,real m_e, real gama);

//Equation 6: P.A variation with time
real alpha_rk(real pper_tmp,real alpha_tmp,real eta_tmp,real m_e, real kz, real w_h, real w, real dwh_ds, real gama, real wtau_sq);

real aeq_rk(real ppar_tmp, real pper_tmp, real alpha_tmp, real eta_tmp, real aeq_tmp, real kappa, real gama, real Bw_out);
