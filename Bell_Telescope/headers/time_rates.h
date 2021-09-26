#pragma once
#include "common.h"
#include "constants.h"
#include <cmath>


//Time rates functions prototypes.

//Equation 1: parallel(to B) speed uz=dz/dt, variation with time
real z_rk(real ppar_tmp, real gama);			

//Equation 2: parallel momentum dpz/dt, variation with time
real p_par_rk(real pper_tmp,real eta_tmp, real w_h, real dwh_ds, real gama);
//Equation 2': Overloaded for interaction. 2 last arguments are also passed here
real p_par_rk(real pper_tmp,real eta_tmp, real w_h, real dwh_ds, real gama, real kz, real wtau_sq);

//Equation 3: differentiate perpendicular momentum, variation with time
real p_per_rk(real ppar_tmp,real pper_tmp,real eta_tmp, real w_h, real dwh_ds, real gama);
//Equation 3': Overlaoded for interaction
real p_per_rk(real ppar_tmp,real pper_tmp,real eta_tmp, real w_h, real dwh_ds, real gama, real w1, real w2, real R1, real R2, real beta);

//Equation 4: eta: angle between BwR and u_per, variation with time
real eta_rk(real ppar_tmp, real w_h, real gama);
//Equation 4':  Overlaoded for interaction
real eta_rk(real ppar_tmp, real w_h, real gama,real kz);

//Equation 5: lamda variation with time
real lamda_rk(real ppar_tmp,real lamda_tmp, real gama);

//Equation 6: P.A variation with time
real alpha_rk(real pper_tmp,real w_h, real dwh_ds, real gama);
//Equation 6': Overloaded for interaction 4 last arguments are also passed here
real alpha_rk(real pper_tmp, real w_h, real dwh_ds, real gama, real alpha_tmp, real eta_tmp, real kz, real wtau_sq);

//Equation 7: Only for interaction. Equatorial P.A variation with time
real aeq_rk(real ppar_tmp, real pper_tmp, real alpha_tmp, real eta_tmp, real aeq_tmp, real kappa, real gama, real Bw_out);
