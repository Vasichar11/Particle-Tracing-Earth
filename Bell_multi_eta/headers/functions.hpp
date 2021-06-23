//Stix tuple, needs wps,wc and returns SDPRL
std::tuple<real, real, real, real, real> stix_parameters(real w, real wce, real wcO, real wcH, real wcHe, real wpse, real wpsO, real wpsH, real wpsHe){ 
//B0mag magnitude of the ambient magnetic field
    real R=1-(wpse/(w*(w+wce)))-(wpsH/(w*(w+wcH)))-(wpsHe/(w*(w+wcHe)))-(wpsO/(w*(w+wcO)));
    real L=1-(wpse/(w*(w-wce)))-(wpsH/(w*(w-wcH)))-(wpsHe/(w*(w-wcHe)))-(wpsO/(w*(w-wcO)));
    real P=1-(wpse/pow(w,2))-(wpsH/pow(w,2))-(wpsHe/pow(w,2))-(wpsO/pow(w,2));
    real S = (R+L)/2;
    real D = (R-L)/2 ;
return std::make_tuple(S, D, P, R, L); //try std::vector for performance
}


//Magnetic dipole field
real Bmag_dipole(real lamda) {                              //Components of B [Walt,1994]
    //----calculate the dipole magnetic field strength      //B = sqrt(Br^2 + Bl^2) = B0*(Re/r^3)*slat_term
    //lamda geomagnetic latitude                            //where r = Re*L*sqrt(clat)
    real slat=sin(lamda);                                   //substituting in Bmag
    real clat=cos(lamda);
    real slat_term = sqrt(1 + 3*slat*slat);             
    real Bmag=((Constants::B0)/(pow((Constants::L_shell),3)))*slat_term/pow(clat,6);
return Bmag;  
}

//Solve dispersion relation
std::tuple<real, real, real, real> dispersion(real S, real P, real R, real L, real D, real th, real w_wave){
    //----- solve dispersion relation (find refractive index and k vector)
    //th wave normal angle
    real A=S*sin(th)*sin(th)+P*cos(th)*cos(th);
    real B=R*L*sin(th)*sin(th)+S*P*(1+cos(th)*cos(th));
    real C=P*R*L; real mu_sq;                       //using Lorentz & Maxwell equations with Stix parameters
    if(B>0)                                         //produced dispersion relation  A*(mu^4) − B*(mu^2) + C = 0
      { mu_sq=(B-sqrt(B*B-4*A*C))/(2*A); }
    if(B<0)
      { mu_sq=(2*C)/(B+sqrt(B*B-4*A*C)); }
    real mu=sqrt(mu_sq);
    real kappa=mu*w_wave/Constants::c;
    real kx=kappa*sin(th);
    real kz=kappa*cos(th);
return std::make_tuple(mu, kappa, kx, kz); //try std::vector for performance
}

//Estimate dwh_ds
real dwh_dsf(real w_h, real lamda){         //[Tao et al, 2012]
    real slat = sin(lamda);
    real clat = cos(lamda);
    real slat_term = sqrt(1 + 3*slat*slat);
    real dwh_ds = (3.0*w_h/(Constants::L_shell*Constants::Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat));
    // dwh_ds = (3.0*w_h/(L_shell*Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term))
return dwh_ds;
}
    
//Estimate Bell parameters //In[10]:                //Bell[1984]
void Bell_params(const real m_e, const real q_e, const real ppar, const real pper, const real Bxw, const real Byw, const real Exw, const real Eyw, const real Ezw, const real kz, const real kx, const real w_h, real &gama, real &w1, real &w2, real &wtau_sq, real &R1, real &R2, real &beta) {
    real p_mag=sqrt(ppar*ppar+pper*pper);
    gama=sqrt((p_mag*p_mag*(Constants::c)*(Constants::c))+(m_e*m_e*(Constants::c)*(Constants::c)*(Constants::c)*(Constants::c)))/(m_e*(Constants::c)*(Constants::c));
    w1=(q_e/(2*m_e))*(Bxw+Byw); 
    w2=(q_e/(2*m_e))*(Bxw-Byw);
    real wtau0_sq=(w1*kz*pper)/(gama*m_e);  //Borntik thesis 2.25d
    beta=(kx*pper)/(m_e*gama*w_h);          //Borntik thesis 2.25a
    real a1=w2/w1;                          //Borntik thesis 2.25f
    real a2=(q_e*Ezw)/(w1*pper);            //Borntik thesis 2.25g
    R1=(Exw+Eyw)/(Bxw+Byw);                 //Borntik thesis 2.25h
    R2=(Exw-Eyw)/(Bxw-Byw);                 //Borntik thesis 2.25h

    wtau_sq = (pow((-1),(Constants::m_res-1)) * wtau0_sq * 
            ( jn( ((Constants::m_res)-1), beta ) - 
                a1*jn( ((Constants::m_res)+1) , beta ) +
                gama*a2*jn( Constants::m_res , beta ) )) ;
//Ji first kind Bessel function
//integer order i=(m_res-1)||(m_res+1)||mres
//cylinder function because we have integer order, solving problem in cylindrical coordinate system
//Borntik thesis 2.25c
//compile with c++17, bessel is defined in <cmath>
//which bessel should i use?
return; 
}

//Estimate resonant velocity
void vres_f(const real w, const real kz, const real w_h, const real alpha, real &v_para_res, real &E_res) {
    real t1 = w*w*kz*kz;
    real t2 = pow((Constants::m_res)*w_h, 2)-(w*w);
    real t3 = kz*kz + pow((Constants::m_res*w_h),2)/(pow((Constants::c)*cos(alpha),2));
    real direction;                 //positive direction counter-streaming particle, negative is co-streaming
    if((Constants::m_res)==0){
        direction=-1.0*copysign(1,kz);}                             //unsigned m_res                     
    else{                                                           
        direction = copysign(1,kz)*copysign(1,Constants::m_res);}   //signed m_res
    v_para_res = ( direction*sqrt(t1 + (t2*t3)) - w*kz) / t3;
    real v_tot_res = v_para_res / cos(alpha);
    E_res = Constants::e_el*(1.0/sqrt( 1-(v_tot_res*v_tot_res/(Constants::c*Constants::c)) ) - 1 );  //not sure about that... taken from Bortnik code
return;
}