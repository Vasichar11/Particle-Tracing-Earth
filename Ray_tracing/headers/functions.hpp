
//Magnetic dipole field #In[19]:
real Bmag_dipole(real lamda)                                      //Components of B [Walt,1994]
{   //----calculate the dipole magnetic field strength            //B = sqrt(Br^2 + Bl^2) = B0*(Re/r^3)*slat_term
    //lamda geomagnetic latitude                                  //Where r = Re*L*sqrt(clat)
    real slat=sin(lamda);                                         //Substituting in Bmag
    real clat=cos(lamda);                   
    real slat_term = sqrt(1 + 3*slat*slat);             
    real Bmag=((Constants::B0)/(pow((Constants::L_shell),3)))*slat_term/pow(clat,6);
    return Bmag;  
}


//Estimate resonant velocity #In[21]:
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
    real E_res = Constants::e_el*(1.0/sqrt( 1-(v_tot_res*v_tot_res/(Constants::c*Constants::c)) ) - 1 ); 
    return v_para_res; //Return by reference, global variables, temp.
}