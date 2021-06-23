#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <tuple>
#include <vector> 
#include <string>
#include <chrono>
#include <fstream> 
#include <iomanip>  // std::setprecision()
//#include <boost/math/special_functions/bessel.hpp>	//std::cyl_bessel_j didn't work well. throws domain error exception in some cases.
//Same directory headers							    //Still boost's bessel gives 14 decimal accuracy relatively py results.
#include "headers/constants.hpp"
#include "headers/functions.hpp"
#include "headers/struct_Particles.h"   //Include libraries(with correct order).
#include "headers/struct_Species.hpp"
#include "headers/struct_Waves.hpp"
//include "/headers/outfile.hpp"  
#include "headers/time_rates.hpp"

//HDF5 I/O
#include <highfive/H5File.hpp>          
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace h5 = HighFive;

int main()
{
	//Objects for each specie.
	Species electron,oxygen,hydrogen,helium;	
	electron.mass    = Constants::m_e;  //electron mass in kg
	electron.charge  = Constants::q_e;  //electron charge in C
	electron.n_factor= 1;	 			//n_factor is needed for density 
	oxygen.mass      = Constants::m_O;  		
	oxygen.charge    = Constants::q_i; 	
	oxygen.n_factor  = 0.006;			
	hydrogen.mass    = Constants::m_H; 		
	hydrogen.charge	 = Constants::q_i; 	
	hydrogen.n_factor= 0.94;
	helium.mass	     = Constants::m_He; 			
	helium.charge    = Constants::q_i; 	
	helium.n_factor  = 0.054;

	//Object for wave.
	Waves wave(Constants::Nsteps + 1);       
	wave.By = Constants::By_wave;  	//Byw=1pT, assign exactly like python
	wave.f = Constants::f_wave; 			  	//wave frequency in Hz
	wave.w = 2*M_PI*wave.f;   	//wave angular frequency 
	wave.th_d = Constants::theta0_deg;
	wave.th = Constants::theta0;//initial wave normal angle.
	
	

	//Find pitch angle at lamda0
	real Beq0=Bmag_dipole(0);
	real Blam0=Bmag_dipole(Constants::lamda0);
	real salpha0=sin(Constants::aeq0)*sqrt(Blam0/Beq0);
	real alpha0=asin(salpha0);
	real w_h_0=(std::abs(electron.charge)*Blam0)/electron.mass;

	//Find momentum from energy
	real Ejoule0=1.602176487E-16*Constants::Ekev0;
	real gama0=(Ejoule0/(electron.mass*pow(Constants::c,2))) + 1;
	real speed0=sqrt( 1 - (1/pow(gama0,2)) ) * Constants::c;
	real zeta0   =  0;
	real upar0=speed0*cos(alpha0);
	real uper0=speed0*sin(alpha0);
	real pper0=gama0*electron.mass*uper0;
	real ppar0=gama0*electron.mass*upar0;



	//Object for particles.
	std::cout<<"\nEta distribution in degrees"<<"\n|From"<<" Mid "<<"To |";
	std::cout<<"\n| "<<Constants::start_d << "  " << Constants::mid_d << " " << Constants:: end_d <<"|\n";

	Particles single(Constants::lamda0, zeta0, upar0, uper0, ppar0, pper0, alpha0, alpha0 ,Constants::aeq0, Constants::aeq0, 0);
	
	std::vector<Particles> eql_dstr(Constants::population, single);
						//equally distributed
	
	//Initialize gyrophases.
	for(int p=0; p<Constants::population; p++)
	{   
	    eql_dstr[p].set_eta((Constants::start_d + p*Constants::lin_step_d)*Constants::D2R);  //RADS
		std::cout<<"\nParticle"<< p << " at " <<  Constants::start_d + p*Constants::lin_step_d ;
	}
	
	//Temp variables that will be assigned in the structure
	real new_lamda , new_zeta, new_uper , new_upar, new_ppar, new_pper, new_alpha, new_alpha2, new_aeq, new_aeq2, new_eta, new_time;

	 


	//Initializations of arrays(some are initialized inside structs)
	std::vector<real> kappa_out;
	std::vector<real> kx_out;
	std::vector<real> kz_out; 
	std::vector<real> L_stix;
	std::vector<real> S_stix;
	std::vector<real> D_stix;
	std::vector<real> P_stix;
	std::vector<real> R_stix;
	std::vector<real> B_earth_out;
	std::vector<real> wh_out;
	std::vector<real> dwh_dt_out;
	std::vector<real> mu_out;
	std::vector<real> mu_adiabatic;
	std::vector<real> vresz_o;
	std::vector<real> Eres_o;
	std::vector<real> E_kin;
	std::vector<real> gama_out;
	std::vector<real> deta_dt;
	std::vector<real> time_sim;
	kappa_out.resize((Constants::Nsteps) + 1);
	kx_out.resize((Constants::Nsteps) + 1);
	kz_out.resize((Constants::Nsteps) + 1);
	L_stix.resize((Constants::Nsteps) + 1);
	S_stix.resize((Constants::Nsteps) + 1);
	D_stix.resize((Constants::Nsteps) + 1);
	P_stix.resize((Constants::Nsteps) + 1);
	R_stix.resize((Constants::Nsteps) + 1);
	B_earth_out.resize((Constants::Nsteps) + 1);
	wh_out.resize((Constants::Nsteps) + 1);
	dwh_dt_out.resize((Constants::Nsteps) + 1);
	mu_out.resize((Constants::Nsteps) + 1);
	mu_adiabatic.resize((Constants::Nsteps) + 1);
	vresz_o.resize((Constants::Nsteps) + 1);
	Eres_o.resize((Constants::Nsteps) + 1);
	E_kin.resize((Constants::Nsteps) + 1);
	gama_out.resize((Constants::Nsteps) + 1);
	deta_dt.resize((Constants::Nsteps) + 1);
	time_sim.resize((Constants::Nsteps) + 1);

	real B0mag, ns_e, wc_e, wps_e, ns_O, wc_O, wps_O ,ns_H, wc_H, wps_H, ns_He, wc_He, wps_He;
	real Bmag,Beq,B_lam,aa,salphaeq,salphaeq2,u_mag;
	real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4;
	real gama,w1,w2,R1,R2,beta,wtau_sq;


	std::cout.precision(10);	//Output 10 decimal precise
	std::cout<<std::scientific;	//For e notation representation
	
	auto rk_start = std::chrono::high_resolution_clock::now();
	

	//Runge kutta
	for(int p=0; p<Constants::population; p++)     //Loop for all particle gyrophases
	{

		int i=0;

		while(i<Constants::Nsteps) 
		{

			Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i));
			B_earth_out.at(i)=Bmag; 
			//Densities
			ns_e = electron.density(eql_dstr[p].lamda.at(i)); ns_O = oxygen.density(eql_dstr[p].lamda.at(i));  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)); ns_He = helium.density(eql_dstr[p].lamda.at(i));
			//Cyclotron frequencies
			wc_e = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);  wc_He = helium.wc(Bmag); wh_out.at(i)=wc_e;
			//Plasma frequencies squared
			wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
			//Stix
			auto stix = stix_parameters(wave.w, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
			S_stix.at(i)=std::get<0>(stix); D_stix.at(i)=std::get<1>(stix); P_stix.at(i)=std::get<2>(stix); R_stix.at(i)=std::get<3>(stix); L_stix.at(i)=std::get<4>(stix);
			//Dispersion
			auto disp = dispersion(S_stix.at(i),P_stix.at(i),R_stix.at(i),L_stix.at(i),D_stix.at(i),Constants::theta0,wave.w);
			mu_out.at(i)=std::get<0>(disp); kappa_out.at(i)=std::get<1>(disp); kx_out.at(i)=std::get<2>(disp); kz_out.at(i)=std::get<3>(disp);
			//Fields
			wave.whistlers(i,mu_out.at(i),P_stix.at(i),D_stix.at(i),S_stix.at(i));
			//Total fields
			wave.Ew.at(i)=sqrt(wave.Exw.at(i)*wave.Exw.at(i) + wave.Eyw.at(i)*wave.Eyw.at(i) + wave.Ezw.at(i)*wave.Ezw.at(i));
			wave.Bw.at(i)=sqrt(wave.Bxw.at(i)*wave.Bxw.at(i) + wave.Byw.at(i)*wave.Byw.at(i) + wave.Bzw.at(i)*wave.Bzw.at(i));
			//dwh_dt
			dwh_dt_out.at(i)=dwh_dsf(wh_out.at(i),eql_dstr[p].lamda.at(i));
			//Bell_params
			Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i),eql_dstr[p].pper.at(i),wave.Bxw.at(i),wave.Byw.at(i),wave.Exw.at(i),wave.Eyw.at(i),wave.Ezw.at(i),kz_out.at(i),kx_out.at(i),wh_out.at(i),gama,w1,w2,wtau_sq,R1,R2,beta);
			//vresonant
			vres_f(wave.w,kz_out.at(i),wh_out.at(i),eql_dstr[p].alpha.at(i),vresz_o.at(i),Eres_o.at(i)); //vres?
			//E kin
			E_kin.at(i)=((gama-1)*electron.mass*Constants::c*Constants::c)*6.2415e15;
			
			//Check
			//std::cout<<"\n" << "lamda " << eql_dstr[p].lamda.at(i)<< "\naeq " << eql_dstr[p].aeq2.at(i) << "\n" ;
			//std::cout<<"\n" << "Bxw " << wave.Bxw.at(i) << "\nByw "<<wave.Byw.at(i)<< "\nBzw "<< wave.Bzw.at(i)<< "\nBw " << wave.Bw.at(i) << "\nExw " << wave.Exw.at(i)<< "\nEyw " << wave.Eyw.at(i) << "\nEzw " << wave.Ezw.at(i) << "\nEw " << wave.Ew.at(i) << "\n";
			//std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He ;
			//std::cout<<"\n" << "L " << L_stix.at(i) << "\nS " <<S_stix.at(i)<< "\nD " << D_stix.at(i) << "\nP " << P_stix.at(i) << "\nR " <<R_stix.at(i) << "\n";
			//std::cout<<"\n" << "mu " << mu_out.at(i) << "\nkappa " << kappa_out.at(i)<< "\nkx " << kx_out.at(i) << "\nkz " <<kz_out.at(i) << "\n";
			//std::cout<<"\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "\n";
			//std::cout<<"\n" << "Eres " << Eres_o.at(i)<< "\nvresz " << vresz_o.at(i) << "\nwtau_sq " << wtau_sq << "\nE_kin " << E_kin.at(i) << "\n";
			//std::cout <<"\n" << "Ekin " << E_kin.at(i);
			//std::cout <<"\n" << "mu_adiabatic" << mu_adiabatic.at(i);
		
			k1=z_rk(eql_dstr[p].ppar.at(i),electron.mass,gama);
			l1=p_par_rk(eql_dstr[p].pper.at(i),eql_dstr[p].eta.at(i),electron.mass,kz_out.at(i),wh_out.at(i),dwh_dt_out.at(i),gama, wtau_sq);
			m1=p_per_rk(eql_dstr[p].ppar.at(i),eql_dstr[p].pper.at(i),eql_dstr[p].eta.at(i),electron.mass,wh_out.at(i),dwh_dt_out.at(i),gama,w1,w2,R1,R2,beta);
			n1=eta_rk(electron.mass,wave.w,kz_out.at(i),eql_dstr[p].ppar.at(i),wh_out.at(i),gama);
			o1=lamda_rk(eql_dstr[p].ppar.at(i),eql_dstr[p].lamda.at(i),electron.mass,gama);
			p1=alpha_rk(eql_dstr[p].pper.at(i), eql_dstr[p].alpha.at(i),eql_dstr[p].eta.at(i), electron.mass,  kz_out.at(i),  wh_out.at(i), wave.w ,dwh_dt_out.at(i),gama, wtau_sq);
			//std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn " << n1<< "\no1 " << o1 << "\np1 " << p1 << "\n";	




			Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);
			ns_e = electron.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1); ns_O = oxygen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1); ns_He = helium.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);
			wc_e = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); wh_out.at(i)=wc_e;
			wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
			stix = stix_parameters(wave.w, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
			S_stix.at(i)=std::get<0>(stix); D_stix.at(i)=std::get<1>(stix); P_stix.at(i)=std::get<2>(stix); R_stix.at(i)=std::get<3>(stix); L_stix.at(i)=std::get<4>(stix);
			disp = dispersion(S_stix.at(i),P_stix.at(i),R_stix.at(i),L_stix.at(i),D_stix.at(i),(Constants::theta0),wave.w);
			mu_out.at(i)=std::get<0>(disp); kappa_out.at(i)=std::get<1>(disp); kx_out.at(i)=std::get<2>(disp); kz_out.at(i)=std::get<3>(disp);
			wave.whistlers(i,mu_out.at(i),P_stix.at(i),D_stix.at(i),S_stix.at(i));
			dwh_dt_out.at(i)=dwh_dsf(wh_out.at(i),eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);
			Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,wave.Bxw.at(i),wave.Byw.at(i),wave.Exw.at(i),wave.Eyw.at(i),wave.Ezw.at(i),kz_out.at(i),kx_out.at(i),wh_out.at(i),gama,w1,w2,wtau_sq,R1,R2,beta);

			k2=z_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,electron.mass,gama);  
			l2=p_par_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1,electron.mass,kz_out.at(i),wh_out.at(i),dwh_dt_out.at(i),gama,wtau_sq);
			m2=p_per_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1,electron.mass,wh_out.at(i),dwh_dt_out.at(i),gama,w1,w2,R1,R2,beta);
			n2=eta_rk(electron.mass,wave.w,kz_out.at(i),eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,wh_out.at(i),gama);
			o2=lamda_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1,electron.mass,gama);
			p2=alpha_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1, eql_dstr[p].alpha.at(i)+0.5*(Constants::h)*p1,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1, electron.mass,  kz_out.at(i),  wh_out.at(i), wave.w ,dwh_dt_out.at(i),gama, wtau_sq);
			//std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\n";
			


			Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);
			ns_e = electron.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2); ns_O = oxygen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2); ns_He = helium.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);
			wc_e = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); wh_out.at(i)=wc_e;
			wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
			stix = stix_parameters(wave.w, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
			S_stix.at(i)=std::get<0>(stix); D_stix.at(i)=std::get<1>(stix); P_stix.at(i)=std::get<2>(stix); R_stix.at(i)=std::get<3>(stix); L_stix.at(i)=std::get<4>(stix);
			disp = dispersion(S_stix.at(i),P_stix.at(i),R_stix.at(i),L_stix.at(i),D_stix.at(i),(Constants::theta0),wave.w);
			mu_out.at(i)=std::get<0>(disp); kappa_out.at(i)=std::get<1>(disp); kx_out.at(i)=std::get<2>(disp); kz_out.at(i)=std::get<3>(disp);
			wave.whistlers(i,mu_out.at(i),P_stix.at(i),D_stix.at(i),S_stix.at(i));
			dwh_dt_out.at(i)=dwh_dsf(wh_out.at(i),eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);
			Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,wave.Bxw.at(i),wave.Byw.at(i),wave.Exw.at(i),wave.Eyw.at(i),wave.Ezw.at(i),kz_out.at(i),kx_out.at(i),wh_out.at(i),gama,w1,w2,wtau_sq,R1,R2,beta);
			
			k3=z_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,electron.mass,gama);  
			l3=p_par_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2,electron.mass,kz_out.at(i),wh_out.at(i),dwh_dt_out.at(i),gama, wtau_sq);
			m3=p_per_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2,electron.mass,wh_out.at(i),dwh_dt_out.at(i),gama,w1,w2,R1,R2,beta);
			n3=eta_rk(electron.mass,wave.w,kz_out.at(i),eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,wh_out.at(i),gama);
			o3=lamda_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2,electron.mass,gama);
			p3=alpha_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2, eql_dstr[p].alpha.at(i)+0.5*(Constants::h)*p2,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2, electron.mass,  kz_out.at(i),  wh_out.at(i), wave.w ,dwh_dt_out.at(i),gama, wtau_sq);
			//std::cout<<"\n" << "k3 " << k3 << "\nl1 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 << "\n";
			


			Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i)+(Constants::h)*o3);
			ns_e = electron.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3); ns_O = oxygen.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3);  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3); ns_He = helium.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3);
			wc_e = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); wh_out.at(i)=wc_e;
			wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
			stix = stix_parameters(wave.w, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
			S_stix.at(i)=std::get<0>(stix); D_stix.at(i)=std::get<1>(stix); P_stix.at(i)=std::get<2>(stix); R_stix.at(i)=std::get<3>(stix); L_stix.at(i)=std::get<4>(stix);
			disp = dispersion(S_stix.at(i),P_stix.at(i),R_stix.at(i),L_stix.at(i),D_stix.at(i),(Constants::theta0),wave.w);
			mu_out.at(i)=std::get<0>(disp); kappa_out.at(i)=std::get<1>(disp); kx_out.at(i)=std::get<2>(disp); kz_out.at(i)=std::get<3>(disp);
			wave.whistlers(i,mu_out.at(i),P_stix.at(i),D_stix.at(i),S_stix.at(i));
			dwh_dt_out.at(i)=dwh_dsf(wh_out.at(i),eql_dstr[p].lamda.at(i)+(Constants::h)*o3);
			Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i)+(Constants::h)*l2,eql_dstr[p].pper.at(i)+(Constants::h)*m3,wave.Bxw.at(i),wave.Byw.at(i),wave.Exw.at(i),wave.Eyw.at(i),wave.Ezw.at(i),kz_out.at(i),kx_out.at(i),wh_out.at(i),gama,w1,w2,wtau_sq,R1,R2,beta);
																						//    l2 here?						
			
			k4=z_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,electron.mass,gama);  
			l4=p_par_rk(eql_dstr[p].pper.at(i)+(Constants::h)*m3,eql_dstr[p].eta.at(i)+(Constants::h)*n3,electron.mass,kz_out.at(i),wh_out.at(i),dwh_dt_out.at(i),gama,wtau_sq);
			m4=p_per_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,eql_dstr[p].pper.at(i)+(Constants::h)*m3,eql_dstr[p].eta.at(i)+(Constants::h)*n3,electron.mass,wh_out.at(i),dwh_dt_out.at(i),gama,w1,w2,R1,R2,beta);
			n4=eta_rk(electron.mass,wave.w,kz_out.at(i),eql_dstr[p].ppar.at(i)+(Constants::h)*l3,wh_out.at(i),gama);
			o4=lamda_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,eql_dstr[p].lamda.at(i)+(Constants::h)*o3,electron.mass,gama);
			p4=alpha_rk(eql_dstr[p].pper.at(i)+(Constants::h)*m3, eql_dstr[p].alpha.at(i)+(Constants::h)*p3,eql_dstr[p].eta.at(i)+(Constants::h)*n3, electron.mass,  kz_out.at(i),  wh_out.at(i), wave.w ,dwh_dt_out.at(i),gama,wtau_sq);
			//std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\n";
			
			new_zeta  = eql_dstr[p].zeta.at(i)+((Constants::h)/6)*(k1+2*k2+2*k3+k4); 
			new_ppar  = eql_dstr[p].ppar.at(i)+((Constants::h)/6)*(l1+2*l2+2*l3+l4);
			new_pper  = eql_dstr[p].pper.at(i)+((Constants::h)/6)*(m1+2*m2+2*m3+m4);
			new_eta   = eql_dstr[p].eta.at(i)+((Constants::h)/6)*(n1+2*n2+2*n3+n4);
			new_lamda = eql_dstr[p].lamda.at(i)+((Constants::h)/6)*(o1+2*o2+2*o3+o4);
			new_alpha = eql_dstr[p].alpha.at(i)+((Constants::h)/6)*(p1+2*p2+2*p3+p4);
			
			deta_dt.at(i+1)=((Constants::h)/6)*(n1+2*n2+2*n3+n4);
			
			new_alpha2 = atan(new_pper/new_ppar);
			aa=atan(new_pper/new_ppar);
			Beq=Bmag_dipole(0);
			B_lam=Bmag_dipole(new_lamda);
			salphaeq=sin(new_alpha)*sqrt(Beq/B_lam);
			
			new_aeq = asin(salphaeq);
			
			salphaeq2=sin(new_alpha2)*sqrt(Beq/B_lam);
			new_aeq2 = asin(salphaeq2);
			
		
			new_upar = new_ppar/(electron.mass*gama);
			new_uper = new_pper/(electron.mass*gama);
			u_mag=sqrt(new_upar+pow(new_uper,2)); //u_mag is right?
			
			mu_adiabatic.at(i+1)=(new_pper*new_pper)/std::abs(wh_out.at(i));
			
			i++;
			time_sim.at(i) = time_sim.at(i-1)+(Constants::h);
			new_time = time_sim.at(i);
	        eql_dstr[p].update_state(new_lamda , new_zeta, new_uper , new_upar, new_ppar, new_pper, new_alpha, new_alpha2, new_aeq, new_aeq2, new_eta, new_time);


			//if(eql_dstr[p].lamda.at(i)>0) {	//Stop at equator,  comment for Non linear WPI
			//	break;}															//happens at l>0
		}				
	}
auto rk_stop = std::chrono::high_resolution_clock::now();  
auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(rk_stop - rk_start);
real rk_time = duration2.count() * pow(10,-6);
std::cout << "\nExecution time " << std::setprecision(8) << rk_time <<" seconds" ;



//-------------------------------------------- OUTPUT DATA HDF5 ------------------------------------------------------//
std::vector<real> eta0;
std::vector<std::vector<real>> aeq2_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1) );
std::vector<std::vector<real>> lamda_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1) );
//ASSIGN FROM STRUCT INTO 2D VECTORS (TEMP)
for(int p=0; p<Constants::population; p++)
{
    eta0.push_back(eql_dstr[p].eta.at(0));  //Initial eta of particles in rads. 1D vector

    for(int j=0; j<Constants::Nsteps; j++)  //Fill 2D vector from vector of structs.
    {
        aeq2_plot[p][j] =  eql_dstr[p].aeq2[j];  
        lamda_plot[p][j] =  eql_dstr[p].lamda[j];  
    }
}
//Create of hdf5 file.
h5::File file("./h5files/particles.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
//Create datasets.
h5::DataSet dataset_aeq2  =  file.createDataSet <real> ("aeq2 2D",         h5::DataSpace::From(aeq2_plot));
h5::DataSet dataset_lamda =  file.createDataSet <real> ("lamda 2D",         h5::DataSpace::From(lamda_plot));
h5::DataSet dataset_eta   =  file.createDataSet <real> ("initial_etas",     h5::DataSpace::From(eta0));
h5::DataSet dataset_time  =  file.createDataSet <real> ("simulation time",  h5::DataSpace::From(time_sim));

dataset_aeq2.write(aeq2_plot);
dataset_lamda.write(lamda_plot);
dataset_eta.write(eta0);
dataset_time.write(time_sim);
//------------------------------------------ OUTPUT DATA HDF5 : END --------------------------------------------------//


//---------------------------------------------- OUTPUT DATA TXT -----------------------------------------------------//
	//Write data in txt files. They will be imported in python code.
	//std::vector<real> plot_lamda, plot_aeq2;
	//plot_lamda.resize((Constants::Nsteps) + 1);
	//plot_aeq2.resize((Constants::Nsteps) + 1);
	//for(i=0;i<(Constants::Nsteps);i++)
	//{
	//	plot_lamda.at(i) = eql_dstr[p].lamda.at(i)/(Constants::RAD); //Vectors divided by RAD
	//	plot_aeq2.at(i) = eql_dstr[p].aeq2.at(i)/(Constants::RAD);
	//}
	//
	//outfile(plot_lamda, plot_aeq2, "./cc_text/aeq.txt");
	//outfile(plot_lamda, plot_aeq2, eql_dstr[p].time, "./cc_text/aeq_time.txt");
	//outfile(plot_lamda, E_kin, "./cc_text/Ekin.txt");
	//outfile(plot_lamda, eql_dstr[p].upar, eql_dstr[p].uper, vresz_o, "./cc_text/vpar_vper_vresz.txt");
	//outfile(plot_lamda, mu_adiabatic, "./cc_text/mu_adiabatic.txt");
	//outfile(plot_lamda, eql_dstr[p].eta, deta_dt, "./cc_text/eta_detadt.txt");
	//outfile(plot_lamda, B_earth_out, "./cc_text/Bearth.txt");
	//outfile(plot_lamda, wave.Bxw, wave.Byw, wave.Bzw, wave.Bw, "./cc_text/Bwave.txt");
	////outfile(plot_lamda, wave.Phi, "./cc_text/Phi.txt");
	//outfile(plot_lamda, wave.Exw, wave.Eyw, wave.Ezw, wave.Ew, "./cc_text/Ewave.txt");
	//outfile(plot_lamda, mu_out, "./cc_text/Refractive_index.txt");
	//outfile(plot_lamda, S_stix,P_stix,D_stix,R_stix,L_stix, "./cc_text/Stix.txt");
	//outfile(plot_lamda, kappa_out, kx_out,kz_out, "./cc_text/kappa.txt");
	//outfile(plot_lamda, wh_out,dwh_dt_out, "./cc_text/wh_dwhdt.txt");
//--------------------------------------------- OUTPUT DATA TXT: END----------------------------------------------------//
	
	//Execution time measurement

	return 0; 
}


