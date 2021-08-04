#include <cmath>
#include <iostream>
#include <cinttypes>  
#include <vector> 
#include <string>
#include <chrono>	//For Execution time.
#include <iomanip>  //For std::setprecision()

//#include <boost/math/special_functions/bessel.hpp>	//std::cyl_bessel_j didn't work well. throws domain error exception in some cases. Still boost's bessel gives 14 decimal accuracy relative to py results.

//Same directory headers							    
//Preprocessor macro instructions are added in files to obey ODR.

#include "headers/common.h"
#include "headers/functions.h"
#include "headers/struct_Particles.h"   		    	
#include "headers/struct_Species.h"					
#include "headers/time_rates.h"
#include "headers/constants.h"

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


//-------------------------------------------------------------DISTRIBUTION OF PARTICLES----------------------------------------------------------------//
	//Object for particles.
	std::cout<<"\nEta distribution in degrees"<<"\n|From"<<" Mid "<<"To |";
	std::cout<<"\n| "<<Constants::eta_start_d << "  " << Constants::eta_mid_d << " " << Constants::eta_end_d <<"|\n";
	std::cout<<"\nWith aeq distribution in degrees"<<"\n|From"<<" Mid "<<"To |";
	std::cout<<"\n| "<<Constants::aeq_start_d << "  " << Constants::aeq_mid_d << " " << Constants::aeq_end_d <<"|\n";
	std::cout<< "\nWith same lamda0_deg: " << Constants::lamda0_deg << "\n";

	Particles single; //Single particle struct.
	std::vector<Particles> eql_dstr(Constants::population, single);	//Vector of structs for particle distribution.
						//equally distributed
	

	//Particle population is divided to sets of (aeq_dstr) particles with different initial aeq.
	//There are (eta_dstr) sets of particles that have different initial eta. Particles within the same set, have the same initial eta.
	real eta0,aeq0;
	real Beq0 = Bmag_dipole(0);
	real Blam0 = Bmag_dipole(Constants::lamda0);

	try //Exception handling. Check if initial parameters are valid.
	{

		for(int e=0, i=0; e<Constants::population; e+=Constants::aeq_dstr, i++)	//e seperates sets.
		{   																	
			for(int a=0; a<Constants::aeq_dstr; a++)							//a counts for different aeqs.				
			{																	
		    	eta0 = (Constants::eta_start_d + i*Constants::eta_step_d)*Constants::D2R; //RADS
		    	aeq0 = (Constants::aeq_start_d + a*Constants::aeq_step_d)*Constants::D2R; //RADS
	
			    eql_dstr[e+a].set_eta(eta0);  		
			    eql_dstr[e+a].set_aeq(aeq0);  
				eql_dstr[e+a].calculations0(Beq0,Constants::lamda0,0,0,0,aeq0,Constants::Ekev0); //Push-back initials: P.A, speed and momentum.
				
				std::cout<<"\nParticle"<< e+a << " with eta0 " <<  eta0*Constants::R2D << ", aeq0 " << aeq0*Constants::R2D << " and alpha0 "<< eql_dstr[e+a].alpha.at(0)*Constants::R2D << " in degrees.\n";
				
			}	
		}
	

		//Test print: initial states of particles
		
		//for(int particle=0;particle<100;particle++)
		//{
		//
		//	std::cout<<"\n" << "particle" << particle << " " <<eql_dstr[particle].lamda.at(0)*Constants::R2D <<" " << eql_dstr[particle].zeta.at(0) <<" " <<eql_dstr[particle].upar.at(0) <<" " << eql_dstr[particle].uper.at(0) <<" " << eql_dstr[particle].ppar.at(0) <<" " << eql_dstr[particle].pper.at(0) <<" " << eql_dstr[particle].alpha.at(0) <<" " << eql_dstr[particle].alpha2.at(0) <<" " << eql_dstr[particle].aeq.at(0)*Constants::R2D <<" " << eql_dstr[particle].aeq2.at(0)*Constants::R2D<<" " << eql_dstr[particle].eta.at(0)*Constants::R2D << " " <<eql_dstr[particle].M_adiabatic.at(0) <<" " << eql_dstr[particle].time.at(0);
		//}

//-------------------------------------------------------------DISTRIBUTION OF PARTICLES----------------------------------------------------------------//

		//Temp variables for Runge Kutta
		real new_lamda , new_zeta, new_uper , new_upar, new_ppar, new_pper, new_alpha, new_alpha2, new_aeq, new_aeq2, new_eta, new_M_adiabatic, new_time;
		real ns_e, wc_e, wps_e, ns_O, wc_O, wps_O ,ns_H, wc_H, wps_H, ns_He, wc_He, wps_He, w_h;
		real Bmag,Beq,B_lam,salphaeq,salphaeq2,u_mag;
		real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;
		real gama,w1,w2,R1,R2,beta,wtau_sq;
		real S,D,P,R,L,mu,dwh_ds,vresz,Eres,kappa,kx,kz;
		real Bxwc,Bzwc,Bywc,Exwc,Eywc,Ezwc,Ewc,Bwc;
		 
		//Initializations of arrays(Attributes of Particles are vectors inside Particles struct).
		//Create vectors instead of big arrays to allocate in heap.
		//Each row(population) has a vector which has cols(Nsteps+1) number of elements, each element initialized to 0.
		std::vector <std::vector<real>> Exw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> Eyw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> Ezw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> Bxw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> Byw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> Bzw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> Bw_out    (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> Ew_out    (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> kappa_out (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> kx_out	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> kz_out	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> L_stix	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> S_stix	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> D_stix	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> P_stix	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> R_stix	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> wh_out    (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> mu_out	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> vresz_o   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> Eres_o	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> gama_out  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> deta_dt   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>> time_sim  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>>dwh_dt_out (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>>B_earth_out(Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		std::vector <std::vector<real>>Ekin       (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
  		//std::vector <std::vector<real>>Phi_out    (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0);


		std::cout.precision(16);	//Output 16 decimal precise
		std::cout<<std::scientific;	//For e notation representation
		auto rk_start = std::chrono::high_resolution_clock::now();

		//Runge kutta
		for(int p=0; p<Constants::population; p++)     //Loop for all particles
		{
			int i=0;

			while(i<Constants::Nsteps) 
			{

				Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i));
				//Densities
				ns_e = electron.density(eql_dstr[p].lamda.at(i)); ns_O = oxygen.density(eql_dstr[p].lamda.at(i));  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)); ns_He = helium.density(eql_dstr[p].lamda.at(i));
				//Cyclotron frequencies
				wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);  wc_He = helium.wc(Bmag);
				//Plasma frequencies squared
				wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
				//Stix
				auto stix = stix_parameters(Constants::w_wave, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
				S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
				//Dispersion
				auto disp = dispersion(S,P,R,L,D,Constants::theta0,Constants::w_wave);
				mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
				//Fields.
				whistlers(p,i,mu,P,D,S,kz,eql_dstr[p].zeta.at(i),time_sim[p][i], Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
				//Total fields.
	    		Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
	    		Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
				//dwh_dt
				dwh_ds=dwh_dsf(w_h,eql_dstr[p].lamda.at(i));
				//Bell_params
				Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i),eql_dstr[p].pper.at(i),Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
				//vresonant
				vres_f(Constants::w_wave,kz,w_h,eql_dstr[p].alpha.at(i),vresz,Eres); //vres?
				//E kin
				Ekin[p][i]=((gama-1)*electron.mass*Constants::c*Constants::c)*6.2415e15;
				B_earth_out[p][i]=Bmag; 
				vresz_o[p][i]= vresz;
	        	Eres_o[p][i] = Eres;
	        	Exw_out[p][i]= Exwc;
	        	Eyw_out[p][i]= Eywc;
	        	Ezw_out[p][i]= Ezwc;
	        	Bxw_out[p][i]= Bxwc;
	        	Byw_out[p][i]= Bywc;
	        	Bzw_out[p][i]= Bzwc;
	        	S_stix[p][i] = S;
	        	D_stix[p][i] = D;
	        	P_stix[p][i] = P;
	        	R_stix[p][i] = R;
	        	L_stix[p][i] = L;
	        	Ew_out[p][i] = Ewc;
	        	Bw_out[p][i] = Bwc;
	        	kx_out[p][i] = kx;
	        	kz_out[p][i] = kz;
	        	wh_out[p][i] = w_h;
	        	mu_out[p][i] = mu;
	        	kappa_out[p][i]  = kappa;
	        	dwh_dt_out[p][i] = dwh_ds;
	        	gama_out[p][i]   = gama;
	        	//Phi_out[p][i]    = Ph_w;
	
				//Check
				//std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He ;
				//std::cout<<"\n" << "lamda " << eql_dstr[p].lamda.at(i)<< "\naeq " << eql_dstr[p].aeq2.at(i) << "\n" ;
				//std::cout<<"\n" << "Bxw " << Bxwc << "\nByw "<<Bywc<< "\nBzw "<<Bzwc<< "\nBw " << Bw_out[p][i] << "\nExw " << Exwc<< "\nEyw " << Eywc << "\nEzw " << Ezwc << "\nEw " << Ew_out[p][i] << "\n";
				//std::cout<<"\n" << "L " << L << "\nS " <<S<< "\nD " << D << "\nP " << P << "\nR " <<R << "\n";
				//std::cout<<"\n" << "mu " << mu << "\nkappa " << kappa<< "\nkx " << kx << "\nkz " <<kz << "\n";
				//std::cout<<"\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "\n";
				//std::cout<<"\n" << "Eres " << Eres<< "\nvresz " << vresz << "\nwtau_sq " << wtau_sq << "\nE_kin " <<  Ekin[p][i] << "\nmu_adiabatic" << eql_dstr[p].M_adiabatic.at(i);
			
				k1=z_rk(eql_dstr[p].ppar.at(i),electron.mass,gama);
				l1=p_par_rk(eql_dstr[p].pper.at(i),eql_dstr[p].eta.at(i),electron.mass,kz,w_h,dwh_ds,gama, wtau_sq);
				m1=p_per_rk(eql_dstr[p].ppar.at(i),eql_dstr[p].pper.at(i),eql_dstr[p].eta.at(i),electron.mass,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
				n1=eta_rk(electron.mass,Constants::w_wave,kz,eql_dstr[p].ppar.at(i),w_h,gama);
				o1=lamda_rk(eql_dstr[p].ppar.at(i),eql_dstr[p].lamda.at(i),electron.mass,gama);
				p1=alpha_rk(eql_dstr[p].pper.at(i),eql_dstr[p].alpha.at(i),eql_dstr[p].eta.at(i), electron.mass, kz, w_h, Constants::w_wave ,dwh_ds,gama, wtau_sq);
				q1=aeq_rk(eql_dstr[p].ppar.at(i),eql_dstr[p].pper.at(i),eql_dstr[p].alpha.at(i), eql_dstr[p].eta.at(i), eql_dstr[p].aeq.at(i), kappa, gama, Bwc); 
				//std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn " << n1<< "\no1 " << o1 << "\np1 " << p1 << "\n";	
//#############################################################################################################################################################################################################################################################################################################//
				Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);
				ns_e = electron.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1); ns_O = oxygen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1); ns_He = helium.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);
				wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
				wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
				stix = stix_parameters(Constants::w_wave, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
				S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
				disp = dispersion(S,P,R,L,D,Constants::theta0,Constants::w_wave);
				mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
				whistlers(p,i,mu,P,D,S,kz,eql_dstr[p].zeta.at(i)+0.5*Constants::h*k1,time_sim[p][i]+0.5*Constants::h*k1, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
				Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
	    		Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
				dwh_ds=dwh_dsf(w_h,eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);
				Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
	
				k2=z_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,electron.mass,gama);  
				l2=p_par_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1,electron.mass,kz,w_h,dwh_ds,gama,wtau_sq);
				m2=p_per_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1,electron.mass,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
				n2=eta_rk(electron.mass,Constants::w_wave,kz,eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,w_h,gama);
				o2=lamda_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1,electron.mass,gama);
				p2=alpha_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1, eql_dstr[p].alpha.at(i)+0.5*(Constants::h)*p1,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1, electron.mass,  kz,  w_h, Constants::w_wave ,dwh_ds,gama, wtau_sq);
				q2=aeq_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,eql_dstr[p].alpha.at(i)+0.5*(Constants::h)*p1, eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1, eql_dstr[p].aeq.at(i)+0.5*(Constants::h)*q1, kappa, gama, Bwc); 
				//std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\n";
//#############################################################################################################################################################################################################################################################################################################//
				Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);
				ns_e = electron.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2); ns_O = oxygen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2); ns_He = helium.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);
				wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
				wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
				stix = stix_parameters(Constants::w_wave, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
				S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
				disp = dispersion(S,P,R,L,D,(Constants::theta0),Constants::w_wave);
				mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
				whistlers(p,i,mu,P,D,S,kz,eql_dstr[p].zeta.at(i)+0.5*Constants::h*k2,time_sim[p][i]*0.5*Constants::h, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
				Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
	    		Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
				dwh_ds=dwh_dsf(w_h,eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);
				Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
				
				k3=z_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,electron.mass,gama);  
				l3=p_par_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2,electron.mass,kz,w_h,dwh_ds,gama, wtau_sq);
				m3=p_per_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2,electron.mass,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
				n3=eta_rk(electron.mass,Constants::w_wave,kz,eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,w_h,gama);
				o3=lamda_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2,electron.mass,gama);
				p3=alpha_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2, eql_dstr[p].alpha.at(i)+0.5*(Constants::h)*p2,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2, electron.mass,  kz,  w_h, Constants::w_wave ,dwh_ds,gama, wtau_sq);
				q3=aeq_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,eql_dstr[p].alpha.at(i)+0.5*(Constants::h)*p2, eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2, eql_dstr[p].aeq.at(i)+0.5*(Constants::h)*q2, kappa, gama, Bwc); 
				//std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 << "\n";
//#############################################################################################################################################################################################################################################################################################################//
				Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i)+(Constants::h)*o3);
				ns_e = electron.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3); ns_O = oxygen.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3);  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3); ns_He = helium.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3);
				wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag);
				wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
				stix = stix_parameters(Constants::w_wave, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
				S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
				disp = dispersion(S,P,R,L,D,(Constants::theta0),Constants::w_wave);
				mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
				whistlers(p,i,mu,P,D,S,kz,eql_dstr[p].zeta.at(i)+Constants::h*k3,time_sim[p][i]+Constants::h, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
				Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
	    		Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
				dwh_ds=dwh_dsf(w_h,eql_dstr[p].lamda.at(i)+(Constants::h)*o3);
				Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i)+(Constants::h)*l2,eql_dstr[p].pper.at(i)+(Constants::h)*m3,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
																							//    l2 here?								
				k4=z_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,electron.mass,gama);  
				l4=p_par_rk(eql_dstr[p].pper.at(i)+(Constants::h)*m3,eql_dstr[p].eta.at(i)+(Constants::h)*n3,electron.mass,kz,w_h,dwh_ds,gama,wtau_sq);
				m4=p_per_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,eql_dstr[p].pper.at(i)+(Constants::h)*m3,eql_dstr[p].eta.at(i)+(Constants::h)*n3,electron.mass,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
				n4=eta_rk(electron.mass,Constants::w_wave,kz,eql_dstr[p].ppar.at(i)+(Constants::h)*l3,w_h,gama);
				o4=lamda_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,eql_dstr[p].lamda.at(i)+(Constants::h)*o3,electron.mass,gama);
				p4=alpha_rk(eql_dstr[p].pper.at(i)+(Constants::h)*m3, eql_dstr[p].alpha.at(i)+(Constants::h)*p3,eql_dstr[p].eta.at(i)+(Constants::h)*n3, electron.mass,  kz,  w_h, Constants::w_wave ,dwh_ds,gama,wtau_sq);
				q4=aeq_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,eql_dstr[p].pper.at(i)+(Constants::h)*m3,eql_dstr[p].alpha.at(i)+(Constants::h)*p3, eql_dstr[p].eta.at(i)+(Constants::h)*n3, eql_dstr[p].aeq.at(i)+(Constants::h)*q3, kappa, gama, Bwc); 
				//std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\n";
//#############################################################################################################################################################################################################################################################################################################//
		
				new_zeta  = eql_dstr[p].zeta.at(i)+((Constants::h)/6)*(k1+2*k2+2*k3+k4); 
				new_ppar  = eql_dstr[p].ppar.at(i)+((Constants::h)/6)*(l1+2*l2+2*l3+l4);
				new_pper  = eql_dstr[p].pper.at(i)+((Constants::h)/6)*(m1+2*m2+2*m3+m4);
				new_eta   = eql_dstr[p].eta.at(i)+((Constants::h)/6)*(n1+2*n2+2*n3+n4);
				new_lamda = eql_dstr[p].lamda.at(i)+((Constants::h)/6)*(o1+2*o2+2*o3+o4);
				new_alpha = eql_dstr[p].alpha.at(i)+((Constants::h)/6)*(p1+2*p2+2*p3+p4);
				new_aeq   = eql_dstr[p].aeq.at(i)+((Constants::h)/6)*(q1+2*q2+2*q3+q4);
				
				deta_dt[p][i+1]=((Constants::h)/6)*(n1+2*n2+2*n3+n4);
				
	
				
				Beq=Bmag_dipole(0);								//Why call it everytime? Isn't it Beq0 always?
				B_lam=Bmag_dipole(new_lamda);
				
				//salphaeq=sin(new_alpha)*sqrt(Beq/B_lam);		//Not reliable.
				//new_aeq = asin(salphaeq);	
				
				new_alpha2 = atan(new_pper/new_ppar);			//Differs a lot with new_alpha?
				salphaeq2=sin(new_alpha2)*sqrt(Beq/B_lam);
				new_aeq2 = asin(salphaeq2);
				
			
				new_upar = new_ppar/(electron.mass*gama);
				new_uper = new_pper/(electron.mass*gama);
				u_mag=sqrt(new_upar+pow(new_uper,2)); //u_mag is right? //not used
				
				new_M_adiabatic=(new_pper*new_pper)/std::abs(w_h);
				
				i++;
				time_sim[p][i] = time_sim[p][i-1]+(Constants::h);
				new_time = time_sim[p][i];
		        eql_dstr[p].update_state(new_lamda , new_zeta, new_uper , new_upar, new_ppar, new_pper, new_alpha, new_alpha2, new_aeq, new_aeq2, new_eta, new_M_adiabatic, new_time);
	
		        //std::cout<<"\n"<< new_alpha << " " << new_zeta << " " << new_ppar<< " " << new_pper<< " " << new_eta << " " <<new_lamda<< " " <<new_aeq ;
				//if(eql_dstr[p].lamda.at(i)>0) {	//Stop at equator, 
				//	break;}															
			}				
		}
	
		auto rk_stop = std::chrono::high_resolution_clock::now();  
		auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(rk_stop - rk_start);
		real rk_time = duration2.count() * pow(10,-6);
		std::cout << "\nExecution time " << std::setprecision(8) << rk_time <<" seconds" ;
	

//------------------------------------------ CATCH EXCEPTIONS : START ------------------------------------------------//

	}	//try end	

	catch(int exception)	 //Temp exceptions. Are these ok?
	{
		if(exception == 99)
		{
			std::cout<< "\n" << "Caught exception: |sinusoidal|>1. P.A out of domain.\n" ;
			return EXIT_FAILURE;	//Will take care of the objects life cycle(destruction)? 
		}
		else if(exception == 98)	
		{
			std::cout<< "\n" << "Caught exception: aeq0=0.\n" ;
			return EXIT_FAILURE;
		}	
	}
	
	catch(std::bad_alloc& ex) //doesn't work.
	{
		std::cout<< "\n Out of Memory Error.\n";
		return EXIT_FAILURE;
	}
//------------------------------------------ CATCH EXCEPTIONS : END -------------------------------------------------//	

//-------------------------------------------- OUTPUT DATA HDF5 ------------------------------------------------------//
	std::vector<real> etas0,aeqs0;
	std::vector<std::vector<real>> aeq2_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>>  aeq_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>>  eta_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>>lamda_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>> time_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	//ASSIGN FROM STRUCT INTO 2D VECTORS (TEMP)
	for(int p=0; p<Constants::population; p++)
	{
	    etas0.push_back(eql_dstr[p].eta.at(0));  //Initial eta of particles in rads. 1D vector
		aeqs0.push_back(eql_dstr[p].aeq.at(0));	//Initial aeq of particles in rads. 1D vector
	    for(unsigned long int i=0; i<eql_dstr[p].lamda.size(); i++)  //Fill 2D vector from vector of structs.
	    {
	        aeq_plot[p][i]  = eql_dstr[p].aeq.at(i);  
	        eta_plot[p][i]  = eql_dstr[p].eta.at(i);  
	        aeq2_plot[p][i] = eql_dstr[p].aeq2.at(i);  
	        lamda_plot[p][i]= eql_dstr[p].lamda.at(i);  
	    	time_plot[p][i] = eql_dstr[p].time.at(i);
	    }
	}
	
	//Create of hdf5 file.
	h5::File file("particles.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	//Create datasets.
	h5::DataSet dataset_aeq2  =  file.createDataSet <real> ("aeq2 2D",           h5::DataSpace::From(aeq2_plot));
	h5::DataSet dataset_aeq   =  file.createDataSet <real> ("aeq 2D",            h5::DataSpace::From(aeq_plot));
	h5::DataSet dataset_eta   =  file.createDataSet <real> ("eta 2D",            h5::DataSpace::From(eta_plot));
	h5::DataSet dataset_lamda =  file.createDataSet <real> ("lamda 2D",          h5::DataSpace::From(lamda_plot));
	h5::DataSet dataset_etas0 =  file.createDataSet <real> ("initial_etas",      h5::DataSpace::From(etas0));
	h5::DataSet dataset_time  =  file.createDataSet <real> ("simulation time",   h5::DataSpace::From(time_plot));
	h5::DataSet dataset_aeqs0 =  file.createDataSet <real> ("initial_aeqs",	     h5::DataSpace::From(aeqs0));
	
	dataset_aeq2.write(aeq2_plot);
	dataset_aeq.write(aeq_plot);
	dataset_eta.write(eta_plot);
	dataset_lamda.write(lamda_plot);
	dataset_etas0.write(etas0);
	dataset_time.write(time_plot);
	dataset_aeqs0.write(aeqs0);

//------------------------------------------ OUTPUT DATA HDF5 : END --------------------------------------------------//
	


return 0; 

}


