#include "headers/read_csv.h"
#include "headers/interpolate.h"
#include "headers/common.h"
#include "headers/constants.h"
#include <vector>
#include <chrono>
#include <iomanip>  

#include <highfive/H5File.hpp>                                      
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace h5 = HighFive;

int main()

{


    auto start1 = std::chrono::high_resolution_clock::now();

//-------------------------------------------------- READ CSV --------------------------------------------------------// 
    //Using function "rcsv" to read Ray data && store values in a single(!) vector.

    std::vector<std::pair<std::string, std::vector<double>>> ray_tracing = rcsv("L2_freq2500_psi-89_lat_0_damping.csv");

    int column_size = ray_tracing.at(0).second.size() ;       //Size of columns

    //Declaring seperate vectors to store Ray data individually.
    std::vector <real> timef, posx, posy, posz, vprelx, vprely, vprelz, vgrelx, vgrely, vgrelz, nx, ny, nz;
    std::vector <real> Bx, By, Bz, w, Ns1, Ns2, Ns3, Ns4, psi, theta_res, Y, Lf, alt, lat, lon, damp;
    std::vector <real> S_stix, D_stix, P_stix, R_stix, L_stix, wmega_e ;
    std::vector <real> Nspec,qs1,qs2,qs3,qs4,ms1,ms2,ms3,ms4,nus1,nus2,nus3,nus4;


    //Allocation of values inside seperate vectors 
    for(int i=0; i<column_size; i++)                
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
//-------------------------------------------------- READ CSV: DONE --------------------------------------------------//
    auto stop1= std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
    real time1 = duration1.count()*pow(10,-6);
    std::cout << "\nReading CSV. Time elapsed: " << time1 <<" seconds" ;

    auto start2 = std::chrono::high_resolution_clock::now();

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


//---------------------------------------------- INTERPOLATE RAY: DONE -----------------------------------------------//
    auto stop2= std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
    real time2 = duration2.count()*pow(10,-6);
    std::cout << "\nInterpolating Ray. Time elapsed: " << time2 <<" seconds" ;

    auto start3 = std::chrono::high_resolution_clock::now();

//------------------------------------------- OTHER VECTOR CALCULATIONS --------------------------------------------- //

  

    //Declare with known size
    std::vector <real> mu_ray, mu_sq, spsi, cpsi, kx_ray, kz_ray, kappa_ray, X_stix;
    std::vector <real> rho1, rho2, Byw_sq, fac1, Byw, Bxw, Bzw, Exw, Eyw, Ezw, Bw_ray, w1, w2, R1, R2 ;

for(size_t i=0; i<time_new.size(); i++)
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
    //std::cout<<"\n"<<Byw.at(i);
    Bw_ray.push_back( sqrt(Bxw.at(i)*Bxw.at(i) + Byw.at(i)*Bzw.at(i) + Byw.at(i)*Bzw.at(i)) );
    //From Bell parameters
    w1.push_back( (Constants::q_e/(2*Constants::m_e))*(Bxw.at(i)+Byw.at(i)) );
    w2.push_back( (Constants::q_e/(2*Constants::m_e))*(Bxw.at(i)-Byw.at(i)) );
    R1.push_back( (Exw.at(i)+Eyw.at(i))/(Bxw.at(i)+Byw.at(i)) );   
    R2.push_back( (Exw.at(i)-Eyw.at(i))/(Bxw.at(i)-Byw.at(i)) );
}
//--------------------------------------- OTHER VECTOR CALCULATIONS: DONE ------------------------------------------- //
    auto stop3 = std::chrono::high_resolution_clock::now();
    auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);
    real time3 = duration3.count()*pow(10,-6);
    std::cout << "\nCalculating other vectors. Time elapsed: " << time3 <<" seconds" ;

    auto start4 = std::chrono::high_resolution_clock::now();


//------------------------------------------------- WRITE HDF5 FILE ------------------------------------------------- //

    h5::File file_out("h5files/interpolated_ray.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);

    //Vectors from interpolation
    //h5::DataSet dataset_posx_int =      file_out.createDataSet("posx_int", posx_int);
    //h5::DataSet dataset_posy_int =      file_out.createDataSet("posy_int", posy_int);
    //h5::DataSet dataset_posz_int =      file_out.createDataSet("posz_int",posz_int);
    //h5::DataSet dataset_vprelx_int =    file_out.createDataSet("vprelx_int", vprelx_int);
    //h5::DataSet dataset_vprely_int =    file_out.createDataSet("vprely_int", vprely_int);
    //h5::DataSet dataset_vprelz_int =    file_out.createDataSet("vprelz_int", vprelz_int);
    //h5::DataSet dataset_vgrelx_int =    file_out.createDataSet("vgrelx_int", vgrelx_int);
    //h5::DataSet dataset_vgrely_int =    file_out.createDataSet("vgrely_int", vgrely_int);
    //h5::DataSet dataset_vgrelz_int =    file_out.createDataSet("vgrelz_int", vgrelz_int);
    //h5::DataSet dataset_nx_int =        file_out.createDataSet("nx_int", nx_int);
    //h5::DataSet dataset_ny_int =        file_out.createDataSet("ny_int", ny_int);
    //h5::DataSet dataset_nz_int =        file_out.createDataSet("nz_int", nz_int);
    //h5::DataSet dataset_Bx_int =        file_out.createDataSet("Bx_int", Bx_int);
    //h5::DataSet dataset_By_int =        file_out.createDataSet("By_int", By_int);
    //h5::DataSet dataset_Bz_int =        file_out.createDataSet("Bz_int", Bz_int);
    //h5::DataSet dataset_w_int =         file_out.createDataSet("w_int", w_int);
    //h5::DataSet dataset_qs1_int =       file_out.createDataSet("qs1_int", qs1_int);
    //h5::DataSet dataset_qs2_int =       file_out.createDataSet("qs2_int", qs2_int);
    //h5::DataSet dataset_qs3_int =       file_out.createDataSet("qs3_int", qs3_int);
    //h5::DataSet dataset_qs4_int =       file_out.createDataSet("qs4_int", qs4_int);
    //h5::DataSet dataset_ms1_int =       file_out.createDataSet("ms1_int", ms1_int);
    //h5::DataSet dataset_ms2_int =       file_out.createDataSet("ms2_int", ms2_int);
    //h5::DataSet dataset_ms3_int =       file_out.createDataSet("ms3_int", ms3_int);
    //h5::DataSet dataset_ms4_int =       file_out.createDataSet("ms4_int", ms4_int);
    //h5::DataSet dataset_Ns1_int =       file_out.createDataSet("Ns1_int", Ns1_int);
    //h5::DataSet dataset_Ns2_int =       file_out.createDataSet("Ns2_int", Ns2_int);
    //h5::DataSet dataset_Ns3_int =       file_out.createDataSet("Ns3_int", Ns3_int);
    //h5::DataSet dataset_Ns4_int =       file_out.createDataSet("Ns4_int", Ns4_int);
    //h5::DataSet dataset_psi_int =       file_out.createDataSet("psi_int", psi_int);
    //h5::DataSet dataset_theta_res_int = file_out.createDataSet("theta_res_int", theta_res_int);
    //h5::DataSet dataset_Y_int =         file_out.createDataSet("Y_int", Y_int);
    //h5::DataSet dataset_Lf_int =        file_out.createDataSet("Lf_int", Lf_int);
    //h5::DataSet dataset_alt_int =       file_out.createDataSet("alt_int", alt_int);
    h5::DataSet dataset_lat_int =       file_out.createDataSet("lat_int", lat_int);
    //h5::DataSet dataset_lon_int =       file_out.createDataSet("lon_int", lon_int);
    //h5::DataSet dataset_damp_int =      file_out.createDataSet("damp_int", damp_int);
    //h5::DataSet dataset_S_stix_int =    file_out.createDataSet("S_stix_int", S_stix_int);
    //h5::DataSet dataset_D_stix_int =    file_out.createDataSet("D_stix_int", D_stix_int);
    //h5::DataSet dataset_P_stix_int =    file_out.createDataSet("P_stix_int", P_stix_int);
    //h5::DataSet dataset_R_stix_int =    file_out.createDataSet("R_stix_int", R_stix_int);
    //h5::DataSet dataset_L_stix_int =    file_out.createDataSet("L_stix_int", L_stix_int);
    //h5::DataSet dataset_wmega_e_int =   file_out.createDataSet("wmega_e_int", wmega_e_int);
    //h5::DataSet dataset_mu_ray =        file_out.createDataSet("mu_ray", mu_ray );
    //h5::DataSet dataset_mu_sq =         file_out.createDataSet("mu_sq", mu_sq );
    //h5::DataSet dataset_spsi =          file_out.createDataSet("spsi", spsi );
    //h5::DataSet dataset_cpsi =          file_out.createDataSet("cpsi", cpsi );
    h5::DataSet dataset_kx_ray =        file_out.createDataSet("kx_ray", kx_ray );
    h5::DataSet dataset_kz_ray =        file_out.createDataSet("kz_ray", kz_ray );
    h5::DataSet dataset_kappa_ray =     file_out.createDataSet("kappa_ray", kappa_ray );
    //h5::DataSet dataset_X_stix =        file_out.createDataSet("X_stix", X_stix );
    //h5::DataSet dataset_rho1 =          file_out.createDataSet("rho1", rho1 );
    //h5::DataSet dataset_rho2 =          file_out.createDataSet("rho2", rho2 );
    //h5::DataSet dataset_Byw_sq =        file_out.createDataSet("Byw_sq", Byw_sq );
    //h5::DataSet dataset_fac1 =          file_out.createDataSet("fac1", fac1 );
    //h5::DataSet dataset_Byw =           file_out.createDataSet("Byw", Byw );
    //h5::DataSet dataset_Bxw =           file_out.createDataSet("Bxw", Bxw );
    h5::DataSet dataset_Bzw =           file_out.createDataSet("Bzw", Bzw );
    //h5::DataSet dataset_Exw =           file_out.createDataSet("Exw", Exw );
    //h5::DataSet dataset_Eyw =           file_out.createDataSet("Eyw", Eyw );
    h5::DataSet dataset_Ezw =           file_out.createDataSet("Ezw", Ezw );
    h5::DataSet dataset_Bw_ray =        file_out.createDataSet("Bw_ray", Bw_ray );
    h5::DataSet dataset_w1 =            file_out.createDataSet("w1", w1 );
    h5::DataSet dataset_w2 =            file_out.createDataSet("w2", w2 );
    h5::DataSet dataset_R1 =            file_out.createDataSet("R1", R1 );
    h5::DataSet dataset_R2 =            file_out.createDataSet("R2", R2 );
 

//----------------------------------------------- WRITE HDF5 FILE: DONE -------------------------------------------- //

    auto stop4= std::chrono::high_resolution_clock::now();
    auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(stop4 - start4);
    real time4 = duration4.count()*pow(10,-6);
    std::cout << "\nWriting in h5 file. Time elapsed: " << time4 <<" seconds\n" ;

    return 0;

}