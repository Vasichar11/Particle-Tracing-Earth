#include "headers/functions.h"


//Return Factors, aeq of runge kutta and kz vector, only when particle in packet.
void f_packet (real index, real &Fpar, real &Fper, real &Ftheta, real &aeq_rk, real &kz, const real pper_tmp, const real ppar_tmp, const real eta_tmp, const real aeqsu_tmp, const real alpha_tmp, const real gama, const real w_h, const real p_mag, const real kx_ray, const real kz_ray, const real kappa_ray, const real Bw_ray, const real Bzw, const real Ezw, const real w1, const real w2, const real R1, const real R2 )
{
    
    kz       =    kz_ray; 
    
    //Using wave's values where the particle is confined
    
    //Calculate beta 
    real beta = (kx_ray*pper_tmp)/(Constants::m_e*gama*w_h); //No need to return.
    
    //Factor calculations for ppar ,pper and theta
    real EL = R1*((2*Constants::m_e)/Constants::q_e)*w1;
    real ER = R2*((2*Constants::m_e)/Constants::q_e)*w2;
    real BR =    ((2*Constants::m_e)/Constants::q_e)*w1;
    real BL =    ((2*Constants::m_e)/Constants::q_e)*w2; 
    
    Fpar   = -(pow((-1),(Constants::m_res)-1)*(-Constants::q_e*(Ezw*jn( (Constants::m_res), beta)
                   +(pper_tmp/(gama*Constants::m_e))*BR*jn( (Constants::m_res)-1, beta)
                   -(pper_tmp/(gama*Constants::m_e))*BL*jn( (Constants::m_res)+1, beta))));
    
    Fper   =  (-Constants::q_e*(ER*jn( (Constants::m_res)-1, beta)+EL*jn( (Constants::m_res)+1, beta)
                   -(ppar_tmp/(gama*Constants::m_e))*BR*jn( (Constants::m_res)-1, beta)
                   +(ppar_tmp/(gama*Constants::m_e))*BL*jn( (Constants::m_res)+1, beta)));
    
    Ftheta = (-Constants::q_e*(ER*jn( (Constants::m_res-1), beta) - EL*jn( (Constants::m_res)+1, beta) 
                   -(ppar_tmp/(gama*Constants::m_e))*BR*jn( (Constants::m_res)-1, beta) 
                   -(ppar_tmp/(gama*Constants::m_e))*BL*jn( (Constants::m_res)+1, beta)
                   +(pper_tmp/(gama*Constants::m_e))*Bzw*jn( (Constants::m_res), beta)));
    
    //Calculate Equatorial P.A change.
    aeq_rk =((Constants::q_e*Bw_ray)/(pow(p_mag,2)))*(tan(aeqsu_tmp)/tan(alpha_tmp))*(((Constants::w_wave/kappa_ray)-(ppar_tmp/(gama*Constants::m_e)))*ppar_tmp-(pow(pper_tmp,2)/(gama*Constants::m_e)))*sin(eta_tmp);

}

//Returns quantities needed, regardless WPI.
void f_always(real &p_mag, real &gama, real &w_h, real &dwh_ds, const real lamda_tmp, const real ppar_tmp, const real pper_tmp)
{   
    
    p_mag  =  sqrt(ppar_tmp*ppar_tmp + pper_tmp*pper_tmp);
    gama   =  sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
    
    real slat = sin(lamda_tmp);   //Bmag and dwh_ds needs                                
    real clat = cos(lamda_tmp);                   
    real slat_term = sqrt(1 + 3*slat*slat); 
    
    real Bmag = ((Constants::B0)/(pow((Constants::L_shell),3)))*slat_term/(pow(clat,6));
    w_h = ((Constants::q_e*Bmag)/Constants::m_e);
    dwh_ds = (3.0*w_h/(Constants::L_shell*Constants::Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat)); // [Tao et al, 2012]
}

//Magnetic dipole field 
real Bmag_dipole(real lamda)                                      //Components of B [Walt,1994]
{   //----calculate the dipole magnetic field strength            //B = sqrt(Br^2 + Bl^2) = B0*(Re/r^3)*slat_term
    //lamda geomagnetic latitude                                  //Where r = Re*L*sqrt(clat)
    real slat=sin(lamda);                                         //Substituting in Bmag
    real clat=cos(lamda);                   
    real slat_term = sqrt(1 + 3*slat*slat);             
    real Bmag=((Constants::B0)/(pow((Constants::L_shell),3)))*slat_term/pow(clat,6);
    return Bmag;  
}

//Estimate resonant velocity
real vres_f(real w_wave, real kz, real w_h, real alpha) 
{
    real t1 = w_wave*w_wave*kz*kz;
    real t2 = pow((Constants::m_res)*w_h, 2)-(w_wave*w_wave);
    real t3 = kz*kz + pow((Constants::m_res*w_h),2)/(pow((Constants::c)*cos(alpha),2));
    real direction;                                                 //Positive direction when counter-streaming particle
    if((Constants::m_res)==0)                                       //Negative is co-streaming
    { direction=-1.0*copysign(1,kz); }                              //Unsigned m_res                     
    else                                                            
    { direction = copysign(1,kz)*copysign(1,Constants::m_res); }    //Signed m_res
    real v_para_res = ( direction*sqrt(t1 + (t2*t3)) - w_wave*kz) / t3;
    real v_tot_res = v_para_res / cos(alpha);
    real E_res = Constants::e_el*(1.0/sqrt( 1-(v_tot_res*v_tot_res/(Constants::c*Constants::c)) ) - 1 ); //not used
    return v_para_res; //Return by reference, global variables, temp.
}