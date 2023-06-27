#include "read_csv.h"
#include "interpolate.h"
#include "common.h"
#include "constants.h"
#include <vector>
#include <chrono>
#include <iomanip>  
#include <iostream>
#include <H5Cpp.h>

int main()

{


    auto start1 = std::chrono::high_resolution_clock::now();

//-------------------------------------------------- READ CSV --------------------------------------------------------// 
    //Using function "rcsv" to read Ray data && store values in a single(!) vector.

    std::vector<std::pair<std::string, std::vector<real>>> ray_tracing = rcsv("src/data/L2_freq2500_psi-89_lat_0_damping.csv");

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
        
    size_t size = timef.size();
    real last_timestep = timef.at(size - 2); // Second to last element.
    int rounded_last_timestep = static_cast<int>(std::round(last_timestep));

    std::vector<real> time_new;
    time_new.reserve(static_cast<size_t>(rounded_last_timestep / Simulation::h));

    for (real t = 0; t < rounded_last_timestep - Simulation::h; t += Simulation::h) {
        time_new.push_back(t);
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
        spsi.push_back( sin(psi_int.at(i)*Universal::D2R) );
        cpsi.push_back( cos(psi_int.at(i)*Universal::D2R) );
        kx_ray.push_back( (nx_int.at(i)*w_int.at(i)) / Universal::c );
        kz_ray.push_back( (nz_int.at(i)*w_int.at(i)) / Universal::c );
        kappa_ray.push_back( (mu_ray.at(i)*w_int.at(i)) / Universal::c );
        X_stix.push_back( P_stix_int.at(i)/(P_stix_int.at(i)-mu_sq.at(i)*spsi.at(i)*spsi.at(i)) );
        rho1.push_back( ((mu_sq.at(i)-S_stix_int.at(i))*mu_sq.at(i)*spsi.at(i)*cpsi.at(i)) / (D_stix_int.at(i)*(mu_sq.at(i)*spsi.at(i)*spsi.at(i)-P_stix_int.at(i))) );
        rho2.push_back( (mu_sq.at(i)-S_stix_int.at(i)) / D_stix_int.at(i) );

        //Define whistler waves #In[7]:R1
        Byw_sq.push_back( ((2.0*Universal::mu_0/Universal::c)*((Wave::pwr*damp_int.at(i))*X_stix.at(i)*X_stix.at(i)*rho2.at(i)*rho2.at(i)*std::abs(cos(psi_int.at(i)*Universal::D2R)))
                        /sqrt(pow((tan(psi_int.at(i)*Universal::D2R)-rho1.at(i)*rho2.at(i)*X_stix.at(i)),2) + pow((1+rho2.at(i)*rho2.at(i)*X_stix.at(i)),2))) );
        fac1.push_back( (P_stix_int.at(i)-mu_sq.at(i)*pow(sin(psi_int.at(i)*Universal::D2R),2)) );
        Byw.push_back( sqrt(Byw_sq.at(i)) );
        Bxw.push_back( std::abs((-(D_stix_int.at(i)*fac1.at(i))/(P_stix_int.at(i)*(S_stix_int.at(i)-mu_sq.at(i))))*Byw.at(i)) );
        Bzw.push_back( std::abs(((D_stix_int.at(i)*sin(psi_int.at(i)*Universal::D2R)*fac1.at(i))/(P_stix_int.at(i)*cos(psi_int.at(i)*Universal::D2R)*(S_stix_int.at(i)-mu_sq.at(i))))*Byw.at(i)) );
        Exw.push_back( std::abs(((Universal::c*fac1.at(i))/(mu_ray.at(i)*P_stix_int.at(i)*cos(psi_int.at(i)*Universal::D2R))*Byw.at(i))) );
        Eyw.push_back( std::abs(((D_stix_int.at(i)*Universal::c*fac1.at(i))/(mu_ray.at(i)*P_stix_int.at(i)*cos(psi_int.at(i)*Universal::D2R)*(pow(mu_ray.at(i),2)-S_stix_int.at(i))))*Byw.at(i)) );
        Ezw.push_back( std::abs((-(Universal::c*mu_ray.at(i)*sin(psi_int.at(i)))/P_stix_int.at(i))*Byw.at(i)) );
        //std::cout<<"\n"<<Byw.at(i);
        Bw_ray.push_back( sqrt(Bxw.at(i)*Bxw.at(i) + Byw.at(i)*Bzw.at(i) + Byw.at(i)*Bzw.at(i)) );
        //From Bell parameters
        w1.push_back( (Universal::q_e/(2*Universal::m_e))*(Bxw.at(i)+Byw.at(i)) );
        w2.push_back( (Universal::q_e/(2*Universal::m_e))*(Bxw.at(i)-Byw.at(i)) );
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

	std::string file_name = "output/files/interpolated_ray_pwr" + std::to_string(Wave::pwr) + ".h5";
	H5::H5File file(file_name, H5F_ACC_TRUNC);
    hsize_t dim = Distribution::population;
	H5::DataSpace data_space(1,&dim);
	H5::DataType data_type(H5::PredType::NATIVE_DOUBLE);


	// Create the datasets
    // Vectors from interpolation
    //H5::DataSet dataset_posx_int =      file.createDataSet("posx_int", data_type, data_space);
    //H5::DataSet dataset_posy_int =      file.createDataSet("posy_int", data_type, data_space);
    //H5::DataSet dataset_posz_int =      file.createDataSet("posz_int",posz_int);
    //H5::DataSet dataset_vprelx_int =    file.createDataSet("vprelx_int", data_type, data_space);
    //H5::DataSet dataset_vprely_int =    file.createDataSet("vprely_int", data_type, data_space);
    //H5::DataSet dataset_vprelz_int =    file.createDataSet("vprelz_int", data_type, data_space);
    //H5::DataSet dataset_vgrelx_int =    file.createDataSet("vgrelx_int", data_type, data_space);
    //H5::DataSet dataset_vgrely_int =    file.createDataSet("vgrely_int", data_type, data_space);
    //H5::DataSet dataset_vgrelz_int =    file.createDataSet("vgrelz_int", data_type, data_space);
    //H5::DataSet dataset_nx_int =        file.createDataSet("nx_int", data_type, data_space);
    //H5::DataSet dataset_ny_int =        file.createDataSet("ny_int", data_type, data_space);
    //H5::DataSet dataset_nz_int =        file.createDataSet("nz_int", data_type, data_space);
    //H5::DataSet dataset_Bx_int =        file.createDataSet("Bx_int", data_type, data_space);
    //H5::DataSet dataset_By_int =        file.createDataSet("By_int", data_type, data_space);
    //H5::DataSet dataset_Bz_int =        file.createDataSet("Bz_int", data_type, data_space);
    //H5::DataSet dataset_w_int =         file.createDataSet("w_int", data_type, data_space);
    //H5::DataSet dataset_qs1_int =       file.createDataSet("qs1_int", data_type, data_space);
    //H5::DataSet dataset_qs2_int =       file.createDataSet("qs2_int", data_type, data_space);
    //H5::DataSet dataset_qs3_int =       file.createDataSet("qs3_int", data_type, data_space);
    //H5::DataSet dataset_qs4_int =       file.createDataSet("qs4_int", data_type, data_space);
    //H5::DataSet dataset_ms1_int =       file.createDataSet("ms1_int", data_type, data_space);
    //H5::DataSet dataset_ms2_int =       file.createDataSet("ms2_int", data_type, data_space);
    //H5::DataSet dataset_ms3_int =       file.createDataSet("ms3_int", data_type, data_space);
    //H5::DataSet dataset_ms4_int =       file.createDataSet("ms4_int", data_type, data_space);
    //H5::DataSet dataset_Ns1_int =       file.createDataSet("Ns1_int", data_type, data_space);
    //H5::DataSet dataset_Ns2_int =       file.createDataSet("Ns2_int", data_type, data_space);
    //H5::DataSet dataset_Ns3_int =       file.createDataSet("Ns3_int", data_type, data_space);
    //H5::DataSet dataset_Ns4_int =       file.createDataSet("Ns4_int", data_type, data_space);
    //H5::DataSet dataset_psi_int =       file.createDataSet("psi_int", data_type, data_space);
    //H5::DataSet dataset_theta_res_int = file.createDataSet("theta_res_int", data_type, data_space);
    //H5::DataSet dataset_Y_int =         file.createDataSet("Y_int", data_type, data_space);
    //H5::DataSet dataset_Lf_int =        file.createDataSet("Lf_int", data_type, data_space);
    //H5::DataSet dataset_alt_int =       file.createDataSet("alt_int", data_type, data_space);
    H5::DataSet dataset_lat_int =       file.createDataSet("lat_int", data_type, data_space);
    //H5::DataSet dataset_lon_int =       file.createDataSet("lon_int", data_type, data_space);
    //H5::DataSet dataset_damp_int =      file.createDataSet("damp_int", data_type, data_space);
    //H5::DataSet dataset_S_stix_int =    file.createDataSet("S_stix_int", data_type, data_space);
    //H5::DataSet dataset_D_stix_int =    file.createDataSet("D_stix_int", data_type, data_space);
    //H5::DataSet dataset_P_stix_int =    file.createDataSet("P_stix_int", data_type, data_space);
    //H5::DataSet dataset_R_stix_int =    file.createDataSet("R_stix_int", data_type, data_space);
    //H5::DataSet dataset_L_stix_int =    file.createDataSet("L_stix_int", data_type, data_space);
    //H5::DataSet dataset_wmega_e_int =   file.createDataSet("wmega_e_int", data_type, data_space);
    //H5::DataSet dataset_mu_ray =        file.createDataSet("mu_ray", data_type, data_space); 
    //H5::DataSet dataset_mu_sq =         file.createDataSet("mu_sq", data_type, data_space); 
    //H5::DataSet dataset_spsi =          file.createDataSet("spsi", data_type, data_space); 
    //H5::DataSet dataset_cpsi =          file.createDataSet("cpsi", data_type, data_space); 
    H5::DataSet dataset_kx_ray =        file.createDataSet("kx_ray", data_type, data_space); 
    H5::DataSet dataset_kz_ray =        file.createDataSet("kz_ray", data_type, data_space); 
    H5::DataSet dataset_kappa_ray =     file.createDataSet("kappa_ray", data_type, data_space); 
    //H5::DataSet dataset_X_stix =        file.createDataSet("X_stix", data_type, data_space); 
    //H5::DataSet dataset_rho1 =          file.createDataSet("rho1", data_type, data_space); 
    //H5::DataSet dataset_rho2 =          file.createDataSet("rho2", data_type, data_space); 
    //H5::DataSet dataset_Byw_sq =        file.createDataSet("Byw_sq", data_type, data_space); 
    //H5::DataSet dataset_fac1 =          file.createDataSet("fac1", data_type, data_space); 
    //H5::DataSet dataset_Byw =           file.createDataSet("Byw", data_type, data_space); 
    //H5::DataSet dataset_Bxw =           file.createDataSet("Bxw", data_type, data_space); 
    H5::DataSet dataset_Bzw =           file.createDataSet("Bzw", data_type, data_space); 
    //H5::DataSet dataset_Exw =           file.createDataSet("Exw", data_type, data_space); 
    //H5::DataSet dataset_Eyw =           file.createDataSet("Eyw", data_type, data_space); 
    H5::DataSet dataset_Ezw =           file.createDataSet("Ezw", data_type, data_space); 
    H5::DataSet dataset_Bw_ray =        file.createDataSet("Bw_ray", data_type, data_space); 
    H5::DataSet dataset_w1 =            file.createDataSet("w1", data_type, data_space); 
    H5::DataSet dataset_w2 =            file.createDataSet("w2", data_type, data_space); 
    H5::DataSet dataset_R1 =            file.createDataSet("R1", data_type, data_space); 
    H5::DataSet dataset_R2 =            file.createDataSet("R2", data_type, data_space); 
    H5::DataSet dataset_Byw =           file.createDataSet("Byw", data_type, data_space);
    H5::DataSet dataset_time =          file.createDataSet("time", data_type, data_space); 
    

	// Write
    //dataset_posx_int.write(posx_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_posy_int.write(posy_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_posz_int.write(posz_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_vprelx_int.write(vprelx_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_vprely_int.write(vprely_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_vprelz_int.write(vprelz_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_vgrelx_int.write(vgrelx_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_vgrely_int.write(vgrely_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_vgrelz_int.write(vgrelz_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_nx_int.write(nx_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_ny_int.write(ny_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_nz_int.write(nz_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Bx_int.write(Bx_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_By_int.write(By_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Bz_int.write(Bz_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_w_int.write(w_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_qs1_int.write(qs1_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_qs2_int.write(qs2_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_qs3_int.write(qs3_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_qs4_int.write(qs4_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_ms1_int.write(ms1_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_ms2_int.write(ms2_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_ms3_int.write(ms3_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_ms4_int.write(ms4_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Ns1_int.write(Ns1_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Ns2_int.write(Ns2_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Ns3_int.write(Ns3_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Ns4_int.write(Ns4_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_psi_int.write(psi_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_theta_res_int.write(theta_res_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Y_int.write(Y_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Lf_int.write(Lf_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_alt_int.write(alt_int.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_lat_int.write(lat_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_lon_int.write(lon_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_damp_int.write(damp_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_S_stix_int.write(S_stix_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_D_stix_int.write(D_stix_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_P_stix_int.write(P_stix_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_R_stix_int.write(R_stix_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_L_stix_int.write(L_stix_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_wmega_e_int.write(wmega_e_int.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_mu_ray.write(mu_ray.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_mu_sq.write(mu_sq.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_spsi.write(spsi.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_cpsi.write(cpsi.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_kx_ray.write(kx_ray.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_kz_ray.write(kz_ray.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_kappa_ray.write(kappa_ray.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_X_stix.write(X_stix.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_rho1.write(rho1.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_rho2.write(rho2.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Byw_sq.write(Byw_sq.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_fac1.write(fac1.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_Byw.write(Byw.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Bxw.write(Bxw.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_Bzw.write(Bzw.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Exw.write(Exw.data(), H5::PredType::NATIVE_DOUBLE);
    //dataset_Eyw.write(Eyw.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_Ezw.write(Ezw.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_Bw_ray.write(Bw_ray.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_w1.write(w1.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_w2.write(w2.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_R1.write(R1.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_R2.write(R2.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_Byw.write(Byw.data(), H5::PredType::NATIVE_DOUBLE);
    dataset_time.write(time_new.data(), H5::PredType::NATIVE_DOUBLE);
//----------------------------------------------- WRITE HDF5 FILE: DONE -------------------------------------------- //

    auto stop4= std::chrono::high_resolution_clock::now();
    auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(stop4 - start4);
    real time4 = duration4.count()*pow(10,-6);
    std::cout << "\nWriting in h5 file. Time elapsed: " << time4 <<" seconds\n" ;

    return 0;

}