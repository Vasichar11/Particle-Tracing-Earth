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
std::tuple<real, real, real, real> dispersion(real S, real P, real R, real L, real D){
    //----- solve dispersion relation (find refractive index and k vector)
    //th wave normal angle
    real A=S*sin(Constants::theta0)*sin(Constants::theta0)+P*cos(Constants::theta0)*cos(Constants::theta0);
    real B=R*L*sin(Constants::theta0)*sin(Constants::theta0)+S*P*(1+cos(Constants::theta0)*cos(Constants::theta0));
    real C=P*R*L;                                   //using Lorentz & Maxwell equations with Stix parameters
//  if(B>0)                                         //produced dispersion relation  A*(mu^4) − B*(mu^2) + C = 0
//    { mu_sq=(B-sqrt(B*B-4*A*C))/(2*A); }
//  if(B<0)
//    { mu_sq=(2*C)/(B+sqrt(B*B-4*A*C)); }
    real mu_sq=(B-sqrt(B*B-4*A*C))/(2*A);
    real mu=sqrt(mu_sq);
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

