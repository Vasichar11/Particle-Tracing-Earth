#include "headers/time_rates.h"


//Function to calculate all needed slopes for the particular step
//(ADIABATIC MOTION)
void slopes(real &k, real &l, real &m, real &n, real &o, real &p, real ppar_tmp, real pper_tmp, real lamda_tmp, real eta_tmp, real w_h, real dwh_ds, real gama)
{
//parallel(to B) speed uz=dz/dt, variation with time
    k = ppar_tmp/(gama*Constants::m_e);

//parallel momentum dpz/dt variation with time
    l = -(1/(Constants::m_e*gama))*((pper_tmp*pper_tmp)/(2*w_h))*dwh_ds; 

//perpendicular momentum variation with time
    m = ((1/(Constants::m_e*gama))*(pper_tmp*ppar_tmp)/(2*w_h))*dwh_ds;

//eta: angle between BwR and u_per variation with time
    n = (((Constants::m_res)*w_h)/gama)-Constants::w_wave;

//lamda variation with time
    o = ppar_tmp/(gama*Constants::m_e*Constants::L_shell*Constants::Re*sqrt(1+3*sin(lamda_tmp)*sin(lamda_tmp))* cos(lamda_tmp));

//P.A variation with time
    p = ((1/(Constants::m_e*gama))*(pper_tmp/(2*w_h)))*dwh_ds;
}


//Overloaded Function to calculate all needed slopes for the particular step
//(WPI)
void slopes(real &k, real &l, real &m, real &n, real &o, real &p, real &q, real ppar_tmp, real pper_tmp, real lamda_tmp, real eta_tmp, real alpha_tmp, real aeq_tmp, real p_mag, real w_h, real dwh_ds, real gama, real kz, real kappa, real wtau_sq, real w1, real w2, real R1, real R2, real beta, real Bw_out)
{
    k = ppar_tmp/(gama*Constants::m_e);   
    l = ((wtau_sq*Constants::m_e)/kz)*sin(eta_tmp)-(1/(Constants::m_e*gama))*((pper_tmp*pper_tmp)/(2*w_h))*dwh_ds;     
    m = -(w1*((ppar_tmp/gama)+(Constants::m_e*R1))*jn((Constants::m_res)-1,beta)
	    - w2*((ppar_tmp/gama)-(Constants::m_e*R2))*jn((Constants::m_res)+1,beta))*sin(eta_tmp)+((1/(Constants::m_e*gama))*(pper_tmp*ppar_tmp)/(2*w_h))*dwh_ds;
    n = (((Constants::m_res)*w_h)/gama)-Constants::w_wave-kz*(ppar_tmp/(Constants::m_e*gama));    
    o = ppar_tmp/(gama*Constants::m_e*Constants::L_shell*Constants::Re*sqrt(1+3*sin(lamda_tmp)*sin(lamda_tmp))* cos(lamda_tmp));   
    p = (-((Constants::m_e*wtau_sq/(kz*pper_tmp))*(1+((cos(alpha_tmp)*cos(alpha_tmp))/(Constants::m_res*(w_h/Constants::w_wave)-1))))*sin(eta_tmp)+
            ((1/(Constants::m_e*gama))*(pper_tmp/(2*w_h)))*dwh_ds);  
//equatorial P.A variation with time
    q  = ((Constants::q_e*Bw_out)/(pow(p_mag,2)))*(tan(aeq_tmp)/tan(alpha_tmp))*(((Constants::w_wave/kappa)-(ppar_tmp/(gama*Constants::m_e)))*ppar_tmp-(pow(pper_tmp,2)/(gama*Constants::m_e)))*sin(eta_tmp);
}