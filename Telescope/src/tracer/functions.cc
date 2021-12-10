#include "headers/functions.h"



//Stix tuple, needs wps,wc and returns SDPRL
std::tuple<real, real, real, real, real> stix_parameters(real wce, real wcO, real wcH, real wcHe, real wpse, real wpsO, real wpsH, real wpsHe){ 
//B0mag magnitude of the ambient magnetic field
    real R=1-(wpse/(Constants::w_wave*(Constants::w_wave+wce)))-(wpsH/(Constants::w_wave*(Constants::w_wave+wcH)))-(wpsHe/(Constants::w_wave*(Constants::w_wave+wcHe)))-(wpsO/(Constants::w_wave*(Constants::w_wave+wcO)));
    real L=1-(wpse/(Constants::w_wave*(Constants::w_wave-wce)))-(wpsH/(Constants::w_wave*(Constants::w_wave-wcH)))-(wpsHe/(Constants::w_wave*(Constants::w_wave-wcHe)))-(wpsO/(Constants::w_wave*(Constants::w_wave-wcO)));
    real P=1-(wpse/pow(Constants::w_wave,2))-(wpsH/pow(Constants::w_wave,2))-(wpsHe/pow(Constants::w_wave,2))-(wpsO/pow(Constants::w_wave,2));
    real S = (R+L)/2;
    real D = (R-L)/2 ;
    return std::make_tuple(S, D, P, R, L); //try std::vector for performance
}

//Magnetic dipole field 
real Bmag_dipole(real lamda)                                      //Components of B [Walt,1994]
{   //----calculate the dipole magnetic field strength            //B = sqrt(Br^2 + Bl^2) = B0*(Re/r^3)*slat_term
    //lamda geomagnetic latitude                                  //Where r = Re*L*(clat^2)
    real slat=sin(lamda);                                         //Substituting in Bmag
    real clat=cos(lamda);                   
    real slat_term = sqrt(1 + 3*slat*slat);             
    real Bmag=((Constants::B0)/(pow((Constants::L_shell),3)))*slat_term/pow(clat,6);
    return Bmag;  
}

//Solve dispersion relation
std::tuple<real, real, real, real> dispersion(real S, real P, real R, real L, real D){
    //----- solve dispersion relation (find refractive index and k vector)
    //th wave normal angle
    real A=S*sin(Constants::theta0)*sin(Constants::theta0)+P*cos(Constants::theta0)*cos(Constants::theta0);
    real B=R*L*sin(Constants::theta0)*sin(Constants::theta0)+S*P*(1+cos(Constants::theta0)*cos(Constants::theta0));
    real C=P*R*L;                                   //using Lorentz & Maxwell equations with Stix parameters //produced dispersion relation  A*(mu^4) − B*(mu^2) + C = 0
    real mu_sq, mu;
    //std::cout<<"\nB "<< B << " A " << A << " C " << C <<"\nsqrtof "<<(B*B-4*A*C);
    //mu squared can become negative. Probably it has to do with the ion densities. Domain error.
    if(B>=0)  mu_sq=(B-sqrt(B*B-4*A*C))/(2*A);   
    else      mu_sq=(2*C)/(B+sqrt(B*B-4*A*C));      //Conditional forms are only for computational accuracy[Kimura, 1966].
    //mu_sq=(B-sqrt(B*B-4*A*C))/(2*A);
    //std::cout<<"\nmu_sq "<< mu_sq<<" mu "<< mu;
    mu=sqrt(mu_sq);
    real kappa=mu*Constants::w_wave/Constants::c;
    real kx=kappa*sin(Constants::theta0);
    real kz=kappa*cos(Constants::theta0);
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
void Bell_params(const real ppar, const real pper, const real Bxw, const real Byw, const real Exw, const real Eyw, const real Ezw, const real kz, const real kx, const real w_h, real gama, real &w1, real &w2, real &wtau_sq, real &R1, real &R2, real &beta) {
    w1=(Constants::q_e/(2*Constants::m_e))*(Bxw+Byw); 
    w2=(Constants::q_e/(2*Constants::m_e))*(Bxw-Byw);
    real wtau0_sq=(w1*kz*pper)/(gama*Constants::m_e);  //Borntik thesis 2.25d
    beta=(kx*pper)/(Constants::m_e*gama*w_h);          //Borntik thesis 2.25a
    real a1=w2/w1;                          //Borntik thesis 2.25f
    real a2=(Constants::q_e*Ezw)/(w1*pper);            //Borntik thesis 2.25g 
    R1=(Exw+Eyw)/(Bxw+Byw);                 //Borntik thesis 2.25h
    R2=(Exw-Eyw)/(Bxw-Byw);                 //Borntik thesis 2.25h
    wtau_sq = (pow((-1),(Constants::m_res-1)) * wtau0_sq * 
            ( jn( ((Constants::m_res)-1), beta ) - 
                a1*jn( ((Constants::m_res)+1) , beta ) +
                gama*a2*jn( Constants::m_res , beta ) )) ;
//Borntik thesis 2.25c

//Bessel function may introduce some inaccuracy(13th decimal)
//All arguments(gama,beta,a1,a2,wtau0_sq) are exactly the same, jn() is not.  
    
//Ji first Kind Bessel function
//integer Order i=(m_res-1)||(m_res+1)||mres
//cylinder function because we have integer order, solving problem in cylindrical coordinate system
//https://en.cppreference.com/w/cpp/numeric/special_functions/cyl_bessel_j
//compile with c++17, bessel is defined in <cmath>
//which bessel should i use?
   return; 
}

//Estimate resonant velocity
void vres_f(const real kz, const real w_h, const real alpha, real &v_para_res, real &E_res) {
    real t1 = Constants::w_wave*Constants::w_wave*kz*kz;
    real t2 = pow((Constants::m_res)*w_h, 2)-(Constants::w_wave*Constants::w_wave);
    real t3 = kz*kz + pow((Constants::m_res*w_h),2)/(pow((Constants::c)*cos(alpha),2));
    real direction;                 //positive direction -> counter-streaming particle, negative-> co-streaming
    if((Constants::m_res)==0){
        direction=-1.0*copysign(1,kz);}                             //unsigned m_res                     
    else{                                                           
        direction = copysign(1,kz)*copysign(1,Constants::m_res);}   //signed m_res
    v_para_res = ( direction*sqrt(t1 + (t2*t3)) - Constants::w_wave*kz) / t3;
    real v_tot_res = v_para_res / cos(alpha);
    E_res = Constants::e_el*(1.0/sqrt( 1-(v_tot_res*v_tot_res/(Constants::c*Constants::c)) ) - 1 );  //not sure about that... taken from Bortnik code
    return;
}

//Compute field components.([Bell 1984] 𝑧̂ ||𝐵0→ and 𝑥̂  pointing towards higher L-shells)
void whistlers(int64_t p, int64_t i, real mu, real P, real D, real S, real kz, real &Bxw,real &Byw,real &Bzw, real &Exw,real &Eyw,real &Ezw)
{

    real mu_sq = pow(mu,2); //mu*mu inaccurate
    real theta = Constants::theta0;
    real fac1 = P-(mu_sq*pow(sin(theta),2)) ;
    Byw = Constants::By_wave;

    Bxw = std::abs((-(D*fac1)/(P*(S-mu_sq)))*Byw); 
    Bzw = std::abs(((D*sin(theta)*fac1)/(P*cos(theta)*(S-(mu_sq)))*Byw));
    Exw = std::abs((((Constants::c)*fac1)/(mu*P*cos(theta))*Byw));
    Eyw = std::abs(((D*(Constants::c)*fac1)/(mu*P*cos(theta)*(mu_sq-S)))*Byw);
    Ezw = std::abs((-((Constants::c)*mu*sin(theta))/P)*Byw);

    //real Phi  = Constants::w_wave*time-kz*zeta; 
    //Exw = -Exw.at(i)*sin(Phi);  
    //Eyw =  Eyw.at(i)*cos(Phi);
    //Ezw = -Ezw.at(i)*sin(Phi);
    //Bxw =  Bxw.at(i)*cos(Phi);  
    //Byw =  Byw.at(i)*sin(Phi);
    //Bzw = -Bzw.at(i)*cos(Phi);

    return;
}






//For Ray Tracing


//Return Factors, aeq of runge kutta and kz vector, only when particle in packet.
void f_packet (real &Fpar, real &Fper, real &Ftheta, real &aeq_rk, real &kz, const real pper_tmp, const real ppar_tmp, const real eta_tmp, const real aeqsu_tmp, const real alpha_tmp, const real gama, const real w_h, const real p_mag, const real kx_ray, const real kz_ray, const real kappa_ray, const real Bw_ray, const real Bzw, const real Ezw, const real w1, const real w2, const real R1, const real R2 )
{
    //Takes and returns same value. Change that
    kz       =    kz_ray; 
    
    //Using wave's values where the particle is confined...
    

    //std::cout<<"\npper_tmp "<<pper_tmp<<"\nppar_tmp "<<ppar_tmp<<"\neta_tmp "<<eta_tmp<<"\naeqsu_tmp "<<aeqsu_tmp<<"\nalpha_tmp "<<alpha_tmp<<"\ngama "<<gama<<"\nw_h "<<w_h<<"\np_mag "<<p_mag<<"\nkx_ray "<<kx_ray<<"\nkz_ray "<<kz_ray<<"\nkappa_ray"<<kappa_ray<<"\nBw_ray "<<Bw_ray<<"\nEzw "<<Ezw<<"\nw1 "<<w1<<"\nw2 "<<w2<<"\nR1 "<<R1<<"\nR2 "<<R2<<"\n";


    //Calculate beta 
    real beta = (kx_ray*pper_tmp)/(Constants::m_e*gama*w_h); //[Bortnic thesis]. No need to return.
    //std::cout<<"\nbeta "<<beta;     //It's negative-> domain error for std::cyl_bessel_j().

    
    //Factor calculations for ppar ,pper and theta
    real EL = R1*((2*Constants::m_e)/Constants::q_e)*w1;
    real ER = R2*((2*Constants::m_e)/Constants::q_e)*w2;
    real BR = ((2*Constants::m_e)/Constants::q_e)*w1;
    real BL = ((2*Constants::m_e)/Constants::q_e)*w2;
    
    Fpar = -(pow((-1),(Constants::m_res-1))*(-Constants::q_e*(Ezw*jn( (Constants::m_res), beta)+(pper_tmp/(gama*Constants::m_e))*BR*jn( (Constants::m_res-1), beta)
                   -(pper_tmp/(gama*Constants::m_e))*BL*jn( (Constants::m_res+1), beta))));
                   
    Fper = (-Constants::q_e*(ER*jn( (Constants::m_res-1), beta)+EL*jn( (Constants::m_res+1), beta)
                  -(ppar_tmp/(gama*Constants::m_e))*BR*jn( (Constants::m_res-1), beta)
                  +(ppar_tmp/(gama*Constants::m_e))*BL*jn( (Constants::m_res+1), beta)));

    Ftheta = (-Constants::q_e*(ER*jn( (Constants::m_res-1), beta) - EL*jn( (Constants::m_res+1), beta) 
                  -(ppar_tmp/(gama*Constants::m_e))*BR*jn( (Constants::m_res-1), beta) 
                  -(ppar_tmp/(gama*Constants::m_e))*BL*jn( (Constants::m_res+1), beta)
                  +(pper_tmp/(gama*Constants::m_e))*Bzw*jn( (Constants::m_res), beta)));
    
    //Calculate Equatorial P.A change.
    aeq_rk =((Constants::q_e*Bw_ray)/(pow(p_mag,2)))*(tan(aeqsu_tmp)/tan(alpha_tmp))*(((Constants::w_wave/kappa_ray)-(ppar_tmp/(gama*Constants::m_e)))*ppar_tmp-(pow(pper_tmp,2)/(gama*Constants::m_e)))*sin(eta_tmp);
    return;
}

//Returns quantities needed, regardless WPI.
void f_always(real &p_mag, real &gama, real &w_h, real &dwh_ds, const real lamda_tmp, const real ppar_tmp, const real pper_tmp)
{   
    
    p_mag  =  sqrt(ppar_tmp*ppar_tmp + pper_tmp*pper_tmp);
    gama   =  sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
    //same with gama = sqrt( 1+(u/c)^2 )
    real slat = sin(lamda_tmp);   //Bmag and dwh_ds needs                                
    real clat = cos(lamda_tmp);                   
    real slat_term = sqrt(1 + 3*slat*slat); 
    
    real Bmag = ((Constants::B0)/(pow((Constants::L_shell),3)))*slat_term/(pow(clat,6));
    w_h = ((Constants::q_e*Bmag)/Constants::m_e);
    dwh_ds = (3.0*w_h/(Constants::L_shell*Constants::Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat)); // [Tao et al, 2012]
    return;
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
    //real v_tot_res = v_para_res / cos(alpha);
    //real E_res = Constants::e_el*(1.0/sqrt( 1-(v_tot_res*v_tot_res/(Constants::c*Constants::c)) ) - 1 ); //not used
    return v_para_res; //Return by reference, global variables, temp.
}