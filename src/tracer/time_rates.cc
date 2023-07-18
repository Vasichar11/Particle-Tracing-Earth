#include "time_rates.h"


// Function to calculate all needed slopes for the particular step
// (ADIABATIC MOTION)
void slopes(real &k, real &l, real &m, real &o, real &p, real ppar_tmp, real pper_tmp, real latitude_tmp, real w_h, real dwh_ds, real gama)
{
// parallel(to B) speed uz=dz/dt, variation with time
    k = ppar_tmp/(gama*Universal::m_e);

// parallel momentum dpz/dt variation with time
    l = -(1/(Universal::m_e*gama))*((pper_tmp*pper_tmp)/(2*w_h))*dwh_ds; 

// perpendicular momentum variation with time
    m = ((1/(Universal::m_e*gama))*(pper_tmp*ppar_tmp)/(2*w_h))*dwh_ds;

// latitude variation with time
    o = ppar_tmp/(gama*Universal::m_e*Distribution::L_shell*Universal::Re*sqrt(1+3*sin(latitude_tmp)*sin(latitude_tmp))* cos(latitude_tmp));

// P.A variation with time
    p = ((1/(Universal::m_e*gama))*(pper_tmp/(2*w_h)))*dwh_ds;
}


// Overloaded Function to calculate all needed slopes for the particular step
// (WPI using Bell equations)
void slopes(real &k, real &l, real &m, real &n, real &o, real &p, real &q, real ppar_tmp, real pper_tmp, real latitude_tmp, real eta_tmp, real alpha_tmp, real aeq_tmp, real p_mag, real w_h, real dwh_ds, real gama, real kz, real kappa, real wtau_sq, real w1, real w2, real R1, real R2, real beta, real Bw_out)
{
	// Speed Parallel(to the static magnetic field Bo). 
    k = ppar_tmp/(gama*Universal::m_e);   
    
	// Parallel momentum variation with time. Shouldn't it be f(time,par). Where's parallel momentum dependence?     
    l = ((wtau_sq*Universal::m_e)/kz)*sin(eta_tmp)-(1/(Universal::m_e*gama))*((pper_tmp*pper_tmp)/(2*w_h))*dwh_ds;     
    
	// Perpendicular momentum variation with time       
    m = -(w1*((ppar_tmp/gama)+(Universal::m_e*R1))*jn((Wave::m_res)-1,beta)
	    - w2*((ppar_tmp/gama)-(Universal::m_e*R2))*jn((Wave::m_res)+1,beta))*sin(eta_tmp)+((1/(Universal::m_e*gama))*(pper_tmp*ppar_tmp)/(2*w_h))*dwh_ds;
    
	// eta  variation with time #In[15]:     //eta: Angle between BwR and u_per, w_h local electron gyrofreq.
    n = (((Wave::m_res)*w_h)/gama)-Wave::w_wave-kz*(ppar_tmp/(Universal::m_e*gama));    
    
	// Lamda time rate.    
    o = ppar_tmp/(gama*Universal::m_e*Distribution::L_shell*Universal::Re*sqrt(1+3*sin(latitude_tmp)*sin(latitude_tmp))* cos(latitude_tmp));   
    
	// P.A variation time rate.                   
    p = (-((Universal::m_e*wtau_sq/(kz*pper_tmp))*(1+((cos(alpha_tmp)*cos(alpha_tmp))/(Wave::m_res*(w_h/Wave::w_wave)-1))))*sin(eta_tmp)+
            ((1/(Universal::m_e*gama))*(pper_tmp/(2*w_h)))*dwh_ds);  
    
    // Equatorial P.A variation time rate.  
    q  = ((Universal::q_e*Bw_out)/(pow(p_mag,2)))*(tan(aeq_tmp)/tan(alpha_tmp))*(((Wave::w_wave/kappa)-(ppar_tmp/(gama*Universal::m_e)))*ppar_tmp-(pow(pper_tmp,2)/(gama*Universal::m_e)))*sin(eta_tmp);
}

// Overloaded Function to calculate all needed slopes for the particular step
// (WPI Ray Tracing using Li equations)
void slopes(real &k, real &l, real &m, real &n, real &o, real &p, real &q, const real ppar_tmp, const real pper_tmp, const real alpha_tmp, const real latitude_tmp, const real eta_tmp, const real Fpar, const real Fper, const real Ftheta, const real gama, const real w_h, const real dwh_ds, const real kz, real Bw_ray, real p_mag, real aeqsu_tmp, real kappa_ray)
{	

	// Speed Parallel(to the static magnetic field Bo). 
	k = ppar_tmp/(gama*Universal::m_e);
	
	// Parallel momentum variation with time. Shouldn't it be f(time,par). Where's parallel momentum dependence?     
	l = Fpar*sin(eta_tmp)-(1/(Universal::m_e*gama))*((pper_tmp*pper_tmp)/(2*w_h))*dwh_ds;
	
	// Perpendicular momentum variation with time       
	m = Fper*sin(eta_tmp) + ((1/(Universal::m_e*gama))*(pper_tmp*ppar_tmp)/(2*w_h))*dwh_ds;
	
	// eta  variation with time #In[15]:     //eta: Angle between BwR and u_per, w_h local electron gyrofreq.
	n = ((Wave::m_res*Ftheta)/pper_tmp)*cos(eta_tmp)  +  ((Wave::m_res*w_h)/gama)  -  Wave::w_wave  -  kz*(ppar_tmp/(Universal::m_e*gama));
	
	// Lamda time rate.    
	o = ppar_tmp / (gama * Universal::m_e * Distribution::L_shell * Universal::Re * sqrt(1+3*sin(latitude_tmp)*sin(latitude_tmp)) * cos(latitude_tmp));
	
	// P.A variation time rate.                   
	p = (-(Fpar/pper_tmp)*(1+(cos(alpha_tmp)*cos(alpha_tmp))/(Wave::m_res*(w_h/(gama*Wave::w_wave))-1))*sin(eta_tmp) + ((1/(Universal::m_e*gama))*(pper_tmp/(2*w_h)))*dwh_ds);
    
    // Equatorial P.A variation time rate.
    q =((Universal::q_e*Bw_ray)/(pow(p_mag,2)))*(tan(aeqsu_tmp)/tan(alpha_tmp))*(((Wave::w_wave/kappa_ray)-(ppar_tmp/(gama*Universal::m_e)))*ppar_tmp-(pow(pper_tmp,2)/(gama*Universal::m_e)))*sin(eta_tmp);
}