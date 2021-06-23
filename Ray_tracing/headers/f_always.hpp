void f_always(real &p_mag, real &gama, real &w_h, real &dwh_ds, const real lamda_tmp, const real ppar_tmp, const real pper_tmp)
{	
//f_always returns quantities needed, regardless WPI.
	
p_mag  =  sqrt(ppar_tmp*ppar_tmp + pper_tmp*pper_tmp);

gama   =  sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);

real slat = sin(lamda_tmp);   //Bmag and dwh_ds needs                                
real clat = cos(lamda_tmp);                   
real slat_term = sqrt(1 + 3*slat*slat); 


real Bmag = ((Constants::B0)/(pow((Constants::L_shell),3)))*slat_term/(pow(clat,6));
 

w_h = ((Constants::q_e*Bmag)/Constants::m_e);

//Estimate dwh_ds  [Tao et al, 2012]
dwh_ds = (3.0*w_h/(Constants::L_shell*Constants::Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat));

}


