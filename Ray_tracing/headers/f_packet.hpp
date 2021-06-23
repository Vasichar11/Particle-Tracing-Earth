void f_packet (real index, real &Fpar, real &Fper, real &Ftheta, real &aeq_rk, real &kz, const real pper_tmp, const real ppar_tmp, const real eta_tmp, const real aeqsu_tmp, const real alpha_tmp, const real gama, const real w_h, const real p_mag, const real kx_ray, const real kz_ray, const real kappa_ray, const real Bw_ray, const real Bzw, const real Ezw, const real w1, const real w2, const real R1, const real R2 )
{
//Using wave's values where the particle is confined
     
kz       =    kz_ray; //

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

//jn? or jn()? compiler?

}


