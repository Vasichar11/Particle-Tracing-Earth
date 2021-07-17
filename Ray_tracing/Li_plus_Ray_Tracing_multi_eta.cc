//Standard library header files.
#include <cmath>                        
#include <iostream>                 
#include <tuple>
#include <vector>   
#include <cinttypes>                     //int64_t type
#include <chrono>
#include <iomanip>                       //setprecision()
#include <algorithm>                     //min_element, max_element for std vector. Is it slow?
//Same directory functions and structs.
#include "headers/constants.h"
#include "headers/struct_Particles.h"  
#include "headers/readCSV.h" 
#include "headers/interpolate.h" 
#include "headers/time_rates.h"
#include "headers/functions.h"
#include "headers/is_in_packet.h"
#include "headers/common.h"
//HDF5 I/O
#include <highfive/H5File.hpp>                                      
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
namespace h5 = HighFive;





int main()
{
//-------------------------------------------------- READ CSV --------------------------------------------------------// 
    //Using header function "rcsv" to read Ray data && store values in a single(!) vector.
    std::vector<std::pair<std::string, std::vector<real>>> ray_tracing = rcsv("L2_freq2500_psi-89_lat_0_damping.csv");
    
    //Declaring seperate vectors to store Ray data individually.
    std::vector <real> timef, posx, posy, posz, vprelx, vprely, vprelz, vgrelx, vgrely, vgrelz, nx, ny, nz;
    std::vector <real> Bx, By, Bz, w, Ns1, Ns2, Ns3, Ns4, psi, theta_res, Y, Lf, alt, lat, lon, damp;
    std::vector <real> S_stix, D_stix, P_stix, R_stix, L_stix, wmega_e ;
    std::vector <real> Nspec,qs1,qs2,qs3,qs4,ms1,ms2,ms3,ms4,nus1,nus2,nus3,nus4; //Constant vectors (for now??)
    //Declaring vectors that will be defined-calculated afterwards.
    std::vector <real> mu_ray, mu_sq, spsi, cpsi, kx_ray, kz_ray, kappa_ray, X_stix;
    std::vector <real> rho1, rho2, Byw_sq, fac1, Byw, Bxw, Bzw, Exw, Eyw, Ezw, Bw_ray, w1, w2, R1, R2 ;

    int column_size = ray_tracing.at(0).second.size() ;       //Size of columns

    //Allocation of values inside seperate vectors 
    for(int i=0; i<column_size; i++)      //All columns are same size, like column 0.
    {   
        timef.push_back(ray_tracing.at(0).second.at(i));
        posx.push_back(ray_tracing.at(1).second.at(i));       //Position of wave
        posy.push_back(ray_tracing.at(2).second.at(i));       //
        posz.push_back(ray_tracing.at(3).second.at(i));       //
        vprelx.push_back(ray_tracing.at(4).second.at(i));     //Phase velocity of wave
        vprely.push_back(ray_tracing.at(5).second.at(i));     //
        vprelz.push_back(ray_tracing.at(6).second.at(i));     //
        vgrelx.push_back(ray_tracing.at(7).second.at(i));     //Group velocity of wave
        vgrely.push_back(ray_tracing.at(8).second.at(i));     //
        vgrelz.push_back(ray_tracing.at(9).second.at(i));     //
        nx.push_back(ray_tracing.at(10).second.at(i));        //Refractive index of propagation medium
        ny.push_back(ray_tracing.at(11).second.at(i));        //
        nz.push_back(ray_tracing.at(12).second.at(i));        //
        Bx.push_back(ray_tracing.at(13).second.at(i));        //Magnetic field components of wave
        By.push_back(ray_tracing.at(14).second.at(i));        //
        Bz.push_back(ray_tracing.at(15).second.at(i));        //                                                                        
        w.push_back(ray_tracing.at(16).second.at(i));         //(constant) Angular frequency of wave
        Nspec.push_back(ray_tracing.at(17).second.at(i));     //(constant) Number of species, ions are single charged.
        qs1.push_back(ray_tracing.at(18).second.at(i));       //(constant) Charge of electron particle
        qs2.push_back(ray_tracing.at(19).second.at(i));       //(constant) Charge of positive Ions    
        qs3.push_back(ray_tracing.at(20).second.at(i));       //(constant)
        qs4.push_back(ray_tracing.at(21).second.at(i));       //(constant)
        ms1.push_back(ray_tracing.at(22).second.at(i));       //(constant) Mass of electron particle 
        ms2.push_back(ray_tracing.at(23).second.at(i));       //(constant) Mass of hydrogen particle 
        ms3.push_back(ray_tracing.at(24).second.at(i));       //(constant) Mass of helium particle    
        ms4.push_back(ray_tracing.at(25).second.at(i));       //(constant) Mass of oxygen particle    
        Ns1.push_back(ray_tracing.at(26).second.at(i));       //Density of electrons in radiation belt
        Ns2.push_back(ray_tracing.at(27).second.at(i));       //Density of hydrogen
        Ns3.push_back(ray_tracing.at(28).second.at(i));       //
        Ns4.push_back(ray_tracing.at(29).second.at(i));       //
        nus1.push_back(ray_tracing.at(30).second.at(i));      //(zero) Collision frequency of electron(?) particles
        nus2.push_back(ray_tracing.at(31).second.at(i));      //(zero) 
        nus3.push_back(ray_tracing.at(32).second.at(i));      //(zero)
        nus4.push_back(ray_tracing.at(33).second.at(i));      //(zero)
        psi.push_back(ray_tracing.at(34).second.at(i));       //Wave normal vector
        theta_res.push_back(ray_tracing.at(35).second.at(i)); //Resonance angle
        Y.push_back(ray_tracing.at(36).second.at(i));         //psi - theta_res
        Lf.push_back(ray_tracing.at(37).second.at(i));        //L_shell
        alt.push_back(ray_tracing.at(38).second.at(i));       //Altitude
        lat.push_back(ray_tracing.at(39).second.at(i));       //Latitude
        lon.push_back(ray_tracing.at(40).second.at(i));       //Longitude
        damp.push_back(ray_tracing.at(41).second.at(i));      //damping factor?
        S_stix.push_back(ray_tracing.at(42).second.at(i));    //Stix parameters
        D_stix.push_back(ray_tracing.at(43).second.at(i));    //
        P_stix.push_back(ray_tracing.at(44).second.at(i));    //
        R_stix.push_back(ray_tracing.at(45).second.at(i));    //
        L_stix.push_back(ray_tracing.at(46).second.at(i));    //
        wmega_e.push_back(ray_tracing.at(47).second.at(i));   //Gyrofrequency of electrons
    }                                                         //-----------Interpolation in constant vectors?---------//

    //vector mu_ray_0 = sqrt(nx*nx+ny*ny+nz*nz);  //Arrays, not used now //Refractive index
    //vector kz_ray_0 = nz*w/Constants::c ;
    //vector freq=w/(2*M_PI);
//-------------------------------------------------- READ CSV: DONE --------------------------------------------------//


//------------------------------------------------- INTERPOLATE RAY --------------------------------------------------//
    //std::cout << timef.at(0) << " " << timef.back();                  //Csv's first and last time value.
        
    int t_size = timef.size();      
    real last_timestep = timef.at(t_size - 2);                          //Second to last element.
    last_timestep = (int)(last_timestep*100 + 0.5) ;                    //Integer (value*100 + 0.5)
    last_timestep = last_timestep/100 ;                                 //---> Round last_timestep, 2 decimals.
    //std::cout<< "\nlast timestep " << last_timestep ;                 //Simulation's last timestep.
    
    std::vector <real> time_new {0} ;                                   //Initialize to use .back() in first iteration.
    int j=1;
    while(time_new.back()+(Constants::h) < last_timestep - Constants::h) //Arrange, 0 to last_timestep with stepsize h
    {                                //Why last_timestep - h ??
        time_new.push_back(0+j*Constants::h); 
        j++ ;
    }
    auto int_start = std::chrono::high_resolution_clock::now();
    //Initializing Vectors using a function that returns vector by value, leads to Copy Elision.
    //(NRVO, "Named Return Value Optimization") C++11 vector moves instead of getting copied
    std::vector <real> posx_int      =   interpolate(time_new, timef, posx);
    std::vector <real> posy_int      =   interpolate(time_new, timef, posy);
    std::vector <real> posz_int      =   interpolate(time_new, timef, posz);
    std::vector <real> vprelx_int    =   interpolate(time_new, timef, vprelx);
    std::vector <real> vprely_int    =   interpolate(time_new, timef, vprely);
    std::vector <real> vprelz_int    =   interpolate(time_new, timef, vprelz);
    std::vector <real> vgrelx_int    =   interpolate(time_new, timef, vgrelx);
    std::vector <real> vgrely_int    =   interpolate(time_new, timef, vgrely);
    std::vector <real> vgrelz_int    =   interpolate(time_new, timef, vgrelz);
    std::vector <real> nx_int        =   interpolate(time_new, timef, nx);
    std::vector <real> ny_int        =   interpolate(time_new, timef, ny);
    std::vector <real> nz_int        =   interpolate(time_new, timef, nz);
    std::vector <real> Bx_int        =   interpolate(time_new, timef, Bx);
    std::vector <real> By_int        =   interpolate(time_new, timef, By);
    std::vector <real> Bz_int        =   interpolate(time_new, timef, Bz);
    std::vector <real> w_int         =   interpolate(time_new, timef, w);
    std::vector <real> qs1_int       =   interpolate(time_new, timef, qs1);
    std::vector <real> qs2_int       =   interpolate(time_new, timef, qs2);
    std::vector <real> qs3_int       =   interpolate(time_new, timef, qs3);
    std::vector <real> qs4_int       =   interpolate(time_new, timef, qs4);
    std::vector <real> ms1_int       =   interpolate(time_new, timef, ms1);
    std::vector <real> ms2_int       =   interpolate(time_new, timef, ms2);
    std::vector <real> ms3_int       =   interpolate(time_new, timef, ms3);
    std::vector <real> ms4_int       =   interpolate(time_new, timef, ms4);
    std::vector <real> Ns1_int       =   interpolate(time_new, timef, Ns1);
    std::vector <real> Ns2_int       =   interpolate(time_new, timef, Ns2);
    std::vector <real> Ns3_int       =   interpolate(time_new, timef, Ns3);
    std::vector <real> Ns4_int       =   interpolate(time_new, timef, Ns4);
    std::vector <real> psi_int       =   interpolate(time_new, timef, psi);
    std::vector <real> theta_res_int =   interpolate(time_new, timef, theta_res);
    std::vector <real> Y_int         =   interpolate(time_new, timef, Y);
    std::vector <real> Lf_int        =   interpolate(time_new, timef, Lf);
    std::vector <real> alt_int       =   interpolate(time_new, timef, alt);
    std::vector <real> lat_int       =   interpolate(time_new, timef, lat);
    std::vector <real> lon_int       =   interpolate(time_new, timef, lon);
    std::vector <real> damp_int      =   interpolate(time_new, timef, damp);
    std::vector <real> S_stix_int    =   interpolate(time_new, timef, S_stix);
    std::vector <real> D_stix_int    =   interpolate(time_new, timef, D_stix);
    std::vector <real> P_stix_int    =   interpolate(time_new, timef, P_stix);
    std::vector <real> R_stix_int    =   interpolate(time_new, timef, R_stix);
    std::vector <real> L_stix_int    =   interpolate(time_new, timef, L_stix);
    std::vector <real> wmega_e_int   =   interpolate(time_new, timef, wmega_e);

    auto int_stop= std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(int_stop - int_start);
    real int_time = duration1.count()*pow(10,-6);
    std::cout << "\nFinished with ray interpolation" ;
    //std::cout << "\nTime elapsed: " << std::setprecision(8) << int_time <<" seconds" ;
    //Interpolation not fast, have to change algorithm.
//---------------------------------------------- INTERPOLATE RAY: DONE -----------------------------------------------//

std::cout.precision(32);    //Output precision of 16 decimals
std::cout<<std::scientific; //Representation with e notation


//------------------------------------------- OTHER VECTOR CALCULATIONS --------------------------------------------- //
for(uint i=0; i<time_new.size(); i++)
{
    //Calculate more parameters #In[6]:
    mu_ray.push_back( sqrt(nx_int.at(i)*nx_int.at(i)+ny_int.at(i)*ny_int.at(i)+nz_int.at(i)*nz_int.at(i)) );
    mu_sq.push_back(  mu_ray.at(i)*mu_ray.at(i) );
    spsi.push_back( sin(psi_int.at(i)*Constants::D2R) );
    cpsi.push_back( cos(psi_int.at(i)*Constants::D2R) );
    kx_ray.push_back( (nx_int.at(i)*w_int.at(i)) / Constants::c );
    kz_ray.push_back( (nz_int.at(i)*w_int.at(i)) / Constants::c );
    kappa_ray.push_back( (mu_ray.at(i)*w_int.at(i)) / Constants::c );
    X_stix.push_back( P_stix_int.at(i)/(P_stix_int.at(i)-mu_sq.at(i)*spsi.at(i)*spsi.at(i)) );
    rho1.push_back( ((mu_sq.at(i)-S_stix_int.at(i))*mu_sq.at(i)*spsi.at(i)*cpsi.at(i)) / (D_stix_int.at(i)*(mu_sq.at(i)*spsi.at(i)*spsi.at(i)-P_stix_int.at(i))) );
    rho2.push_back( (mu_sq.at(i)-S_stix_int.at(i)) / D_stix_int.at(i) );

    //Define whistler waves #In[7]:R1
    Byw_sq.push_back( ((2.0*Constants::mu_0/Constants::c)*((Constants::pwr*damp_int.at(i))*X_stix.at(i)*X_stix.at(i)*rho2.at(i)*rho2.at(i)*std::abs(cos(psi_int.at(i)*Constants::D2R)))
                    /sqrt(pow((tan(psi_int.at(i)*Constants::D2R)-rho1.at(i)*rho2.at(i)*X_stix.at(i)),2) + pow((1+rho2.at(i)*rho2.at(i)*X_stix.at(i)),2))) );
    fac1.push_back( (P_stix_int.at(i)-mu_sq.at(i)*pow(sin(psi_int.at(i)*Constants::D2R),2)) );
    Byw.push_back( sqrt(Byw_sq.at(i)) );
    Bxw.push_back( std::abs((-(D_stix_int.at(i)*fac1.at(i))/(P_stix_int.at(i)*(S_stix_int.at(i)-mu_sq.at(i))))*Byw.at(i)) );
    Bzw.push_back( std::abs(((D_stix_int.at(i)*sin(psi_int.at(i)*Constants::D2R)*fac1.at(i))/(P_stix_int.at(i)*cos(psi_int.at(i)*Constants::D2R)*(S_stix_int.at(i)-mu_sq.at(i))))*Byw.at(i)) );
    Exw.push_back( std::abs(((Constants::c*fac1.at(i))/(mu_ray.at(i)*P_stix_int.at(i)*cos(psi_int.at(i)*Constants::D2R))*Byw.at(i))) );
    Eyw.push_back( std::abs(((D_stix_int.at(i)*Constants::c*fac1.at(i))/(mu_ray.at(i)*P_stix_int.at(i)*cos(psi_int.at(i)*Constants::D2R)*(pow(mu_ray.at(i),2)-S_stix_int.at(i))))*Byw.at(i)) );
    Ezw.push_back( std::abs((-(Constants::c*mu_ray.at(i)*sin(psi_int.at(i)))/P_stix_int.at(i))*Byw.at(i)) );

    Bw_ray.push_back( sqrt(Bxw.at(i)*Bxw.at(i) + Byw.at(i)*Bzw.at(i) + Byw.at(i)*Bzw.at(i)) );
    //From Bell parameters
    w1.push_back( (Constants::q_e/(2*Constants::m_e))*(Bxw.at(i)+Byw.at(i)) );
    w2.push_back( (Constants::q_e/(2*Constants::m_e))*(Bxw.at(i)-Byw.at(i)) );
    R1.push_back( (Exw.at(i)+Eyw.at(i))/(Bxw.at(i)+Byw.at(i)) );   
    R2.push_back( (Exw.at(i)-Eyw.at(i))/(Bxw.at(i)-Byw.at(i)) );
}
//--------------------------------------- OTHER VECTOR CALCULATIONS: DONE ------------------------------------------- //




//----------------------------------------- PITCH ANGLE 0 && MOMENTUM 0 --------------------------------------------- //
//Find pitch angle at lamda0 #In[22]:
real Beq0    =  Bmag_dipole(0);
real Blam0   =  Bmag_dipole(Constants::lamda0);          //lamda0 is in rads
real salpha0 =  sin(Constants::aeq0)*sqrt(Blam0/Beq0);
real alpha0  =  asin(salpha0);
real w_h_0   =  (std::abs(Constants::q_e)*Blam0) / Constants::m_e;  //not used
//std::cout << "\n" << Beq0 << " " << Blam0 << " " << Constants::lamda0 ;
//std::cout << "\n" << Constants::aeq0*180/M_PI << "  " << salpha0 << "\n " << alpha0 <<"  "<< alpha0*180/M_PI <<"\n " ;


//Find momentum from energy #In[23]:
real Ejoule0 =  1.602176487E-16*Constants::Ekev0;
real gama0   =  (Ejoule0/(Constants::m_e*pow(Constants::c,2))) + 1;
real speed0  =  sqrt( 1 - (1/pow(gama0,2)) ) * Constants::c;
real zeta0   =  0;
real upar0   =  speed0*cos(alpha0);
real uper0   =  speed0*sin(alpha0);
real pper0   =  gama0*Constants::m_e*uper0;
real ppar0   =  gama0*Constants::m_e*upar0;
real mu_adiabatic_0 = (pper0*pper0) / (2*Constants::m_e*Blam0);
//std::cout << "mu_adiabatic_0" << mu_adiabatic_0 ;
//----------------------------------------- PITCH ANGLE 0 && MOMENTUM 0: DONE ----------------------------------------//





//------------------------------------------- DECLARE VECTORS && VARIABLES -------------------------------------------//
//Declaration of vectors
std::vector <std::vector<real>> zeta    (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> time_sim(Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Exw_out(Constants::population,  std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Eyw_out(Constants::population,  std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Ezw_out(Constants::population,  std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Bxw_out(Constants::population,  std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Byw_out(Constants::population,  std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Bzw_out(Constants::population,  std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Ew_out(Constants::population,   std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Bw_out(Constants::population,   std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> mu_out(Constants::population,   std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> mu_ad_li(Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> vresz_o(Constants::population,  std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Eres_o(Constants::population,   std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Ekin_li(Constants::population,  std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> deta_dt(Constants::population,  std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> kappa_out(Constants::population,std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> wh_out(Constants::population,   std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> eps_res(Constants::population,  std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> ind_reson(Constants::population,std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Ftheta_o(Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Fpar_o(Constants::population,   std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> Fper_o(Constants::population,   std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector <std::vector<real>> gama_out(Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
std::vector<std::vector<real>>mu_adiabatic(Constants::population,std::vector<real> (Constants::Nsteps + 1, 0) );


std::vector <real>::iterator ptr, requested_index;    //To find requested index
int index;                                            //To find minimum difference between latitudes

for(int p=0; p<Constants::population; p++)
{
    Ekin_li[p][0]=Ejoule0;
}


//Variables declaration
real min_lat,max_lat;
int index_res;
int i;
real p_mag,gama,w_h,dwh_ds,kz,Bxw_o,Byw_o,Bzw_o,Exw_o,Eyw_o,Ezw_o,Bw_o,Ew_o,vresz;
real Beq,B_lam,salphaeq,u_mag,Fpar,Fper,Ftheta;
real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;
//--------------------------------------- DECLARE VECTORS && VARIABLES: DONE -----------------------------------------//




//--------------------------- OBJECT CREATION FOR PARTICLES WITH DIFFERENT GYROPHASES --------------------------------//

Particles single(Constants::lamda0, zeta0, upar0, uper0, ppar0, pper0, alpha0, Constants::aeq0, Constants::aeq0,mu_adiabatic_0, 0);

std::vector<Particles> eql_dstr(Constants::population, single);

                    //equally distributed

//Initialize gyrophases.
for(int p=0; p<Constants::population; p++)
{   
    eql_dstr[p].set_eta((0+p*Constants::lin_step_d)*Constants::D2R);  //RADS
}

//Temp variables that will be assigned in the structure
real new_lamda , new_zeta, new_uper , new_upar, new_ppar, new_pper, new_alpha, new_aeq, new_aeqsu, new_eta, new_mu_ad_li, new_time;
//----------------------- OBJECT CREATION FOR PARTICLES WITH DIFFERENT GYROPHASES: DONE ------------------------------//    






//---------------------------------------------------- WPI -----------------------------------------------------------//
auto rk_start = std::chrono::high_resolution_clock::now();

for(int p=0; p<Constants::population; p++)//Loop for all particle gyrophases
{

    i=0;

    while(i<Constants::Nsteps)           //Nsteps-1? 
    {   
        //Take the data from interpolation and return (pulse_dur) of them, moved by "i" each iteration.
        min_lat=*min_element(lat_int.cbegin() + i, lat_int.cbegin() + Constants::puls_dur + i);  //Minimum lat of wave(in pulse duration).
        max_lat=*max_element(lat_int.cbegin() + i, lat_int.cbegin() + Constants::puls_dur + i);  //Max lat of wave.
        //std::cout<<"\n" << eql_dstr[p].lamda.at(i)*Constants::R2D << " " << min_lat << "" << max_lat << " "<< lat_int.at(i) << " " << lat_int.back(); 
         


        f_always(p_mag, gama, w_h, dwh_ds, eql_dstr[p].lamda.at(i), eql_dstr[p].ppar.at(i), eql_dstr[p].pper.at(i));
        Bxw_o = 0;
        Byw_o = 0;
        Bzw_o = 0;
        Exw_o = 0;
        Eyw_o = 0;
        Ezw_o = 0;
        Bw_o  = 0;
        Ew_o  = 0;
        kz    = 0;
        Ftheta= 0; 
        Fpar  = 0;
        Fper  = 0;
        vresz = 0;
        eps_res[p][i] = 1;
        index_res = 0;
        q1 = 0;      
        //is_in_packet returns "if and where" WPI happens.  
        index = is_in_packet(min_lat, max_lat, eql_dstr[p].lamda.at(i), i, lat_int); 
        if(index!=0) { //Do only if there's WPI
            f_packet(index, Fpar, Fper, Ftheta, q1, kz, eql_dstr[p].pper.at(i), eql_dstr[p].ppar.at(i), eql_dstr[p].eta.at(i), eql_dstr[p].aeqsu.at(i), eql_dstr[p].alpha.at(i), gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
            Bxw_o=Bxw[index];// Fields of Whistler wave when WPI is present.
            Byw_o=Byw[index];//
            Bzw_o=Bzw[index];//
            Exw_o=Exw[index];//
            Eyw_o=Eyw[index];//
            Ezw_o=Ezw[index];//
            Bw_o=sqrt(Bxw_o*Bxw_o+Byw_o*Byw_o+Bzw_o*Bzw_o);//
            Ew_o=sqrt(Exw_o*Exw_o+Eyw_o*Eyw_o+Ezw_o*Ezw_o);//
            vresz = vres_f(Constants::w_wave,kz,w_h,eql_dstr[p].alpha.at(i));//
            index_res=1;//
            eps_res[p][i]=std::abs((eql_dstr[p].upar.at(i)-vresz)/vresz); //Maybe one function for eps_res instead of vres_f?
        }
        vresz_o[p][i]=vresz;
        ind_reson[p][i]=index_res;
        Ftheta_o[p][i]=Ftheta;
        Fpar_o[p][i]=Fpar;
        Fper_o[p][i]=Fper;
        Exw_out[p][i]=Exw_o; 
        Eyw_out[p][i]=Eyw_o;
        Ezw_out[p][i]=Ezw_o;
        Bxw_out[p][i]=Bxw_o; 
        Byw_out[p][i]=Byw_o;
        Bzw_out[p][i]=Bzw_o;
        Bw_out[p][i]=Bw_o;
        Ew_out[p][i]=Ew_o;
        gama_out[p][i]=gama;
        slopes(k1,  l1,  m1,  n1,  o1,  p1, eql_dstr[p].ppar.at(i), eql_dstr[p].pper.at(i), eql_dstr[p].alpha.at(i), eql_dstr[p].lamda.at(i), eql_dstr[p].eta.at(i), Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz);
        //Check1
        //std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn1 " << n1<< "\no1 " << o1 << "\np1 " << p1 <<"\nq1 " << q1 <<"\n";
        //std::cout <<"\n"<<"\nBzw"<< Bzw[index] <<"\nw1 "<<w1[index] <<"\nw2 "<< w2[index]  <<"\nR1 "<< R1[index] <<"\nR2 "<< R2[index] <<"\ngama "<< gama <<"\nbeta "<< beta <<"\nppar "<< eql_dstr[p].ppar.at(i) <<"\npper "<< eql_dstr[p].pper.at(i) <<"\nBw_su "<< Bw_su <<"\nFper "<< Fper <<"\nFpar "<< Fpar <<"\nFtheta "<< Ftheta ;  

        f_always(p_mag, gama, w_h, dwh_ds, eql_dstr[p].lamda.at(i)+0.5*Constants::h*o1, eql_dstr[p].ppar.at(i)+0.5*Constants::h*l1, eql_dstr[p].pper.at(i)+0.5*Constants::h*m1);
        kz     =  0;
        Ftheta =  0;
        Fpar   =  0;
        Fper   =  0;
        q2     =  0; 
        index = is_in_packet(min_lat, max_lat, eql_dstr[p].lamda.at(i)+0.5*Constants::h*o1, i, lat_int);
        if(index!=0) { //Do only if there's WPI
            f_packet(index, Fpar, Fper, Ftheta, q2, kz, eql_dstr[p].pper.at(i)+0.5*Constants::h*m1, eql_dstr[p].ppar.at(i)+0.5*Constants::h*l1, eql_dstr[p].eta.at(i)+0.5*Constants::h*n1, eql_dstr[p].aeqsu.at(i)+0.5*Constants::h*q1, eql_dstr[p].alpha.at(i)+0.5*Constants::h*p1, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes( k2, l2, m2, n2, o2, p2, eql_dstr[p].ppar.at(i) + 0.5*Constants::h*l1, eql_dstr[p].pper.at(i) + 0.5*Constants::h*m1, eql_dstr[p].alpha.at(i) + 0.5*Constants::h*p1, eql_dstr[p].lamda.at(i) + 0.5*Constants::h*o1, eql_dstr[p].eta.at(i) + 0.5*Constants::h*n1, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz ); 
        //Check2
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\nq2"<< q2 <<"\n";;
        
        f_always(p_mag, gama, w_h, dwh_ds, eql_dstr[p].lamda.at(i)+0.5*Constants::h*o2, eql_dstr[p].ppar.at(i)+0.5*Constants::h*l2, eql_dstr[p].pper.at(i)+0.5*Constants::h*m2);
        kz     =  0;
        Ftheta =  0;
        Fpar   =  0;
        Fper   =  0;
        q3     =  0; 
        index = is_in_packet(min_lat, max_lat, eql_dstr[p].lamda.at(i)+0.5*Constants::h*o2, i, lat_int);
        if(index!=0) { //Do only if there's WPI
            f_packet(index, Fpar, Fper, Ftheta, q3, kz, eql_dstr[p].pper.at(i)+0.5*Constants::h*m2, eql_dstr[p].ppar.at(i)+0.5*Constants::h*l2, eql_dstr[p].eta.at(i)+0.5*Constants::h*n2, eql_dstr[p].aeqsu.at(i)+0.5*Constants::h*q2, eql_dstr[p].alpha.at(i)+0.5*Constants::h*p2, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes( k3, l3, m3, n3, o3, p3, eql_dstr[p].ppar.at(i) + 0.5*Constants::h*l2, eql_dstr[p].pper.at(i) + 0.5*Constants::h*m2, eql_dstr[p].alpha.at(i) + 0.5*Constants::h*p2, eql_dstr[p].lamda.at(i) + 0.5*Constants::h*o2, eql_dstr[p].eta.at(i) + 0.5*Constants::h*n2, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz );   
        //Check3
        //std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 << "\nq3"<< q3 <<"\n";

        f_always(p_mag, gama, w_h, dwh_ds, eql_dstr[p].lamda.at(i)+Constants::h*o3, eql_dstr[p].ppar.at(i)+Constants::h*l3, eql_dstr[p].pper.at(i)+Constants::h*m3);
        kz     =  0;
        Ftheta =  0;
        Fpar   =  0;
        Fper   =  0;
        q4     =  0; 
        index = is_in_packet(min_lat, max_lat, eql_dstr[p].lamda.at(i)+Constants::h*o3, i, lat_int);
        if(index!=0) { //Do only if there's WPI
            f_packet(index, Fpar, Fper, Ftheta, q4, kz, eql_dstr[p].pper.at(i)+Constants::h*m3, eql_dstr[p].ppar.at(i)+Constants::h*l3, eql_dstr[p].eta.at(i)+Constants::h*n3, eql_dstr[p].aeqsu.at(i)+Constants::h*q3, eql_dstr[p].alpha.at(i)+Constants::h*p3, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes( k4, l4, m4, n4, o4, p4, eql_dstr[p].ppar.at(i) + Constants::h*l3, eql_dstr[p].pper.at(i) + Constants::h*m3, eql_dstr[p].alpha.at(i) + Constants::h*p3, eql_dstr[p].lamda.at(i) + Constants::h*o3, eql_dstr[p].eta.at(i) + Constants::h*n3, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz ); 
        //Check4
        //std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\nq4" << q4 <<"\n";

        //std::cout <<"\np_mag "<<p_mag <<"\ngama "<< gama  <<"\nBmag "<< Bmag <<"\nw_h "<< w_h <<"\ndwh_ds "<< dwh_ds <<"\nbeta "<< beta <<"\nkz "<< kz <<"\nkappa "<< kappa_su <<"\nBw_su "<< Bw_su <<"\nFper "<< Fper <<"\nFpar "<< Fpar <<"\nFtheta "<< Ftheta ;  
        //Fper,Fpar,Ftheta may have 14th decimal diff with python, don't know why, may be Bessel function.(?)

        new_zeta   =  eql_dstr[p].zeta.at(i)   +  (Constants::h/6)*(k1+2*k2+2*k3+k4);
        new_ppar   =  eql_dstr[p].ppar.at(i)   +  (Constants::h/6)*(l1+2*l2+2*l3+l4);
        new_pper   =  eql_dstr[p].pper.at(i)   +  (Constants::h/6)*(m1+2*m2+2*m3+m4);
        new_eta    =  eql_dstr[p].eta.at(i)    +  (Constants::h/6)*(n1+2*n2+2*n3+n4);
        new_lamda  =  eql_dstr[p].lamda.at(i)  +  (Constants::h/6)*(o1+2*o2+2*o3+o4);
        new_alpha  =  eql_dstr[p].alpha.at(i)  +  (Constants::h/6)*(p1+2*p2+2*p3+p4);
        new_aeqsu  =  eql_dstr[p].aeqsu.at(i)  +  (Constants::h/6)*(q1+2*q2+2*q3+q4);
        
        deta_dt[p][i+1] = (Constants::h/6)*(n1+2*n2+2*n3+n4); //is this right ??

        Beq      =  Bmag_dipole(0);
        B_lam    =  Bmag_dipole(new_lamda);
        salphaeq =  sin(new_alpha)*sqrt(Beq/B_lam);
        new_aeq  =  asin(salphaeq);
        
        new_upar = new_ppar / (Constants::m_e*gama);
        new_uper = new_pper / (Constants::m_e*gama);
        u_mag    = sqrt(new_upar+pow(new_uper,2));      //not used
        
        new_mu_ad_li=(new_pper*new_pper)/(2*Constants::m_e*B_lam); 
        
        i++;  
        time_sim[p][i]=time_sim[p][i-1]+Constants::h; 
        new_time = time_sim[p][i];

        eql_dstr[p].update_state(new_lamda , new_zeta, new_uper , new_upar, new_ppar, new_pper, new_alpha, new_aeq, new_aeqsu, new_eta, new_mu_ad_li, new_time);
        //if(eql_dstr[p].lamda.at(i)>0)     //stop at equator
        //    break;


    }

}


auto rk_stop = std::chrono::high_resolution_clock::now();  
auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(rk_stop - rk_start);
real rk_time = duration2.count() * pow(10,-6);
std::cout << "\nFinished with Runge Kutta " ;
std::cout << "\nTime elapsed: " << std::setprecision(8) << rk_time <<" seconds" ;
//------------------------------------------------- WPI: END -------------------------------------------------//


//-------------------------------------------------- OUTPUT DATA HDF5 -----------------------------------------------------//
std::vector<real> eta0;
std::vector<std::vector<real>> aeqsu_plot   (Constants::population, std::vector<real> (Constants::Nsteps + 1) );
std::vector<std::vector<real>> lamda_plot   (Constants::population, std::vector<real> (Constants::Nsteps + 1) );
std::vector<std::vector<real>> mu_ad_li_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1) );
std::vector<std::vector<real>> alpha_plot   (Constants::population, std::vector<real> (Constants::Nsteps + 1) );


//ASSIGN FROM STRUCT INTO 2D VECTORS (TEMP)
for(int p=0; p<Constants::population; p++)
{
    eta0.push_back(eql_dstr[p].eta.at(0));  //Initial eta of particles in rads. 1D vector

    for(int j=0; j<Constants::Nsteps; j++)  //Fill 2D vector from vector of structs.
    {
        aeqsu_plot[p][j]    =  eql_dstr[p].aeqsu[j];  
        lamda_plot[p][j]    =  eql_dstr[p].lamda[j];
        mu_ad_li_plot[p][j] =  eql_dstr[p].mu_ad_li[j];  
        alpha_plot[p][j]    =  eql_dstr[p].alpha[j];
    }
}
//Create of hdf5 file.
h5::File file("particles.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
//Create datasets.
h5::DataSet dataset_aeqsu     =  file.createDataSet <real> ("aeqsu 2D",         h5::DataSpace::From(aeqsu_plot));
h5::DataSet dataset_lamda     =  file.createDataSet <real> ("lamda 2D",         h5::DataSpace::From(lamda_plot));
h5::DataSet dataset_eta       =  file.createDataSet <real> ("initial_etas",     h5::DataSpace::From(eta0));
h5::DataSet dataset_time      =  file.createDataSet <real> ("simulation time",  h5::DataSpace::From(time_sim));
h5::DataSet dataset_mu_ad_li  =  file.createDataSet <real> ("mu_ad_li 2D",      h5::DataSpace::From(mu_ad_li_plot));
h5::DataSet dataset_alpha     =  file.createDataSet <real> ("alpha 2D",         h5::DataSpace::From(alpha_plot));
h5::DataSet dataset_Ftheta    =  file.createDataSet <real> ("Ftheta 2D",        h5::DataSpace::From(Ftheta_o));
h5::DataSet dataset_gama      =  file.createDataSet <real> ("gamma 2D",         h5::DataSpace::From(gama_out));
h5::DataSet dataset_Bxw       =  file.createDataSet <real> ("Bx 2D",            h5::DataSpace::From(Bxw_out));
h5::DataSet dataset_Byw       =  file.createDataSet <real> ("By 2D",            h5::DataSpace::From(Byw_out));
h5::DataSet dataset_Bzw       =  file.createDataSet <real> ("Bz 2D",            h5::DataSpace::From(Bzw_out));


dataset_aeqsu.write(aeqsu_plot);
dataset_lamda.write(lamda_plot);
dataset_eta.write(eta0);
dataset_time.write(time_sim);
dataset_mu_ad_li.write(mu_ad_li_plot);
dataset_alpha.write(alpha_plot);
dataset_Ftheta.write(Ftheta_o);
dataset_gama.write(gama_out);
dataset_Bxw.write(Bxw_out);
dataset_Byw.write(Byw_out);
dataset_Bzw.write(Bzw_out);

//------------------------------------------------- OUTPUT DATA HDF5 : END -------------------------------------------------//

std::cout << "\n\n";
return 0; 

}




