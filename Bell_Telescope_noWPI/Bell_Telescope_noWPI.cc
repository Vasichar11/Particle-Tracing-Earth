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
#include "headers/struct_Telescope.h"  				
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
//-----------------------------------------------------------PARTICLE TELESCOPE PARAMETERS---------------------------------------------------------------//
//Omni-Directional Particle-Telescope.	
	
	//real max_step = 8*pow(10,-5);  // max( |lamda.at(i+1)-lamda.at(i)| )			
	Telescope ODPT(Constants::telescope_lamda, Constants::L_shell); //Satellite situated lat~ 30 deg, L_shell=5. Lamda tolerance = max step of particle. 		
					
	//ODPT.view		    = 360;								//Detector's Angular coverage(degrees).
	//ODPT.sector_range = 15;								//Sector's 	 Angular coverage(degrees).
	//ODPT.sector_num = ( ODPT.view / ODPT.sector_range);   //Number of sectors. Omni->360.	
	//ODPT.tRes    	  = 0.1; 					    		//0.1s time resolution.
//-----------------------------------------------------------PARTICLE TELESCOPE PARAMETERS:END----------------------------------------------------------//

//-------------------------------------------------------------DISTRIBUTION OF PARTICLES----------------------------------------------------------------//
	//Object for particles.
	//std::cout<<"\nEta distribution in degrees"<<"\n|From "<<" To|";
	//std::cout<<"\n| "<<Constants::eta_start_d << "  "<< " " << Constants::eta_end_d <<"|\n";
	std::cout<<"\nWith aeq distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::aeq_start_d << "  "<< " " << Constants::aeq_end_d <<"|\n";
	//std::cout<<"\nWith lamda distribution in degrees"<<"\n|From "<<" To|";
	//std::cout<<"\n| "<<Constants::lamda_start_d << "  "<< " " << Constants::lamda_end_d <<"|\n";

	Particles single; //Single particle struct.
	std::vector<Particles> eql_dstr(Constants::population, single);	//Vector of structs for particle distribution.
						//equally distributed
	
	real eta0,aeq0,lamda0;
	real Beq0 = Bmag_dipole(0);	//Beq isn't always Beq0?

	try //Exception handling. Check if initial parameters are valid.
	{

		for(int e=0, i=0; e<Constants::eta_dstr, i<Constants::population; e++)		//Redundant conditions.
		{   																	
			//eta0 = (Constants::eta_start_d + e*Constants::eta_step_d) * Constants::D2R;			      //RADS
			eta0   = Constants::eta0;	//All same eta0. No WPI.
			
			for(int a=0; a<Constants::aeq_dstr; a++)											
			{	
			    aeq0 = (Constants::aeq_start_d + a*Constants::aeq_step_d) * Constants::D2R;	          //RADS

			    for(int l=0; l<Constants::lamda_dstr; l++, i++)
			    {
			    	//lamda0 = (Constants::lamda_start_d + l*Constants::lamda_step_d)*Constants::D2R;   //RADS  
					lamda0 = Constants::lamda0;	//All same lamda.


				    eql_dstr[i].set_eta(eta0);  		
				    eql_dstr[i].set_aeq(aeq0);  
				    eql_dstr[i].set_lamda(lamda0);  
					eql_dstr[i].calculations0(Beq0,lamda0,0,0,0,aeq0,Constants::Ekev0); //Push-back initials: P.A, speed and momentum.
					
					std::cout<<"\nParticle"<<i<<" with eta0: "<< eta0*Constants::R2D <<", aeq0: "<< aeq0*Constants::R2D <<" and lamda0: "<<lamda0*Constants::R2D<< " in degrees.\n";
				
				}
			}	
		}
		//Test print: initial states of particles
		//for(int particle=0;particle<Constants::population;particle++)
		//{
		//
		//	std::cout<<"\n" << "particle" << particle << " " <<eql_dstr[particle].lamda.at(0)*Constants::R2D <<" " << eql_dstr[particle].zeta.at(0) <<" " <<eql_dstr[particle].upar.at(0) <<" " << eql_dstr[particle].uper.at(0) <<" " << eql_dstr[particle].ppar.at(0) <<" " << eql_dstr[particle].pper.at(0) <<" " << eql_dstr[particle].alpha.at(0) <<" " << eql_dstr[particle].alpha2.at(0) <<" " << eql_dstr[particle].aeq.at(0)*Constants::R2D <<" " << eql_dstr[particle].aeq2.at(0)*Constants::R2D<<" " << eql_dstr[particle].eta.at(0)*Constants::R2D << " " <<eql_dstr[particle].M_adiabatic.at(0) <<" " << eql_dstr[particle].time.at(0);
		//}
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES:END------------------------------------------------------------//

//--------------------------------------------------------------------RUNGE KUTTA-----------------------------------------------------------------------//

		//Temp variables for Runge Kutta
		real new_lamda , new_zeta, new_uper , new_upar, new_ppar, new_pper, new_alpha, new_alpha2, new_aeq, new_aeq2, new_eta, new_M_adiabatic, new_time;
		real ns_e, wc_e, wps_e, ns_O, wc_O, wps_O ,ns_H, wc_H, wps_H, ns_He, wc_He, wps_He, w_h;
		real Bmag,Beq,B_lam,salphaeq,salphaeq2,u_mag;
		real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;
		real gama,w1,w2,R1,R2,beta,wtau_sq;
		real S,D,P,R,L,mu,dwh_ds,vresz,Eres,kappa,kx,kz;
		real Bxwc,Bzwc,Bywc,Exwc,Eywc,Ezwc,Ewc,Bwc;
		real p_mag;
		std::tuple<real, real, real, real, real> stix;
		std::tuple<real, real, real, real> disp;
		real min_diff=100;	
		real max_diff=-100;	

		//Initializations of arrays(Attributes of Particles are vectors inside Particles struct).
		//Create vectors instead of big arrays to allocate in heap.
		//Each row(population) has a vector which has cols(Nsteps+1) number of elements, each element initialized to 0.
		std::vector <std::vector<real>> Exw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
		
		//Temp check size of vectors and if it can be allocated.
		int64_t szv = (Exw_out.size() * Exw_out[0].size() * sizeof(real)) + (Exw_out.size()*sizeof(Exw_out[0])) + sizeof(Exw_out);
		std::cout<< "\nAllocated memory for 2D Vectors is ~: " << szv <<"B * (26+13) = " << szv*39*pow(10,-9) << "GB \n";
		//If these vectors need more than 5GB(for my home system) then abort.
		if(szv*40*pow(10,-9)>5) { throw 97; } 

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

		//std::cout.precision(16);		//Output 16 decimal precise
		//std::cout<<std::scientific;		//For e notation representation
		auto rk_start = std::chrono::high_resolution_clock::now();



		for(int p=0; p<Constants::population; p++)     //Loop for all particles
		{
			int i=0;

			while(i<Constants::Nsteps) 
			{
				Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i));
				//Densities.
				ns_e = electron.density(eql_dstr[p].lamda.at(i)); ns_O = oxygen.density(eql_dstr[p].lamda.at(i));  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)); ns_He = helium.density(eql_dstr[p].lamda.at(i));
				//Cyclotron frequencies.
				wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);  wc_He = helium.wc(Bmag);
				//dwh_dt.
				dwh_ds=dwh_dsf(w_h,eql_dstr[p].lamda.at(i));
				//From bell parameters function(temp for no interaction).
				p_mag = sqrt(eql_dstr[p].ppar.at(i)*eql_dstr[p].ppar.at(i)+eql_dstr[p].pper.at(i)*eql_dstr[p].pper.at(i));
				gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(electron.mass*electron.mass*Constants::c*Constants::c*Constants::c*Constants::c))/(electron.mass*Constants::c*Constants::c);
				
				//Assign no interaction values.
				kz=0;
				wtau_sq=0;
				w1=0;
				w2=0;
				R1=0;
				R2=0;
				beta=0;
				kappa=0;
				Bwc=0;

				B_earth_out[p][i]=Bmag; 
		        wh_out[p][i] = w_h;
		        dwh_dt_out[p][i] = dwh_ds;
		        gama_out[p][i]=gama;
				


				if(Constants::interaction)
				{
					//Plasma frequencies squared.
					wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
					//Stix.
					stix = stix_parameters(Constants::w_wave, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
					S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
					//Dispersion.
					disp = dispersion(S,P,R,L,D,Constants::theta0,Constants::w_wave);
					mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
					//Fields.
					whistlers(p,i,mu,P,D,S,kz,eql_dstr[p].zeta.at(i),time_sim[p][i], Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
					//Total fields.
		    		Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
		    		Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);

					//Bell_params.
					Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i),eql_dstr[p].pper.at(i),Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
					//vresonant.
					vres_f(Constants::w_wave,kz,w_h,eql_dstr[p].alpha.at(i),vresz,Eres); //vres?
					//E kin.
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
		        	mu_out[p][i] = mu;
		        	kappa_out[p][i]  = kappa;
		        	wh_out[p][i] = w_h;
		        	dwh_dt_out[p][i] = dwh_ds;
		        	gama_out[p][i]   = gama;
	        		//Phi_out[p][i]    = Ph_w;
				}
				//Check print parameters
				//std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He << "\nlamda " << eql_dstr[p].lamda.at(i)<< "\naeq " << eql_dstr[p].aeq2.at(i) << "\nBxw " << Bxwc << "\nByw "<<Bywc<< "\nBzw "<<Bzwc<< "\nBw " << Bw_out[p][i] << "\nExw " << Exwc<< "\nEyw " << Eywc << "\nEzw " << Ezwc << "\nEw " << Ew_out[p][i] << "\nL " << L << "\nS " <<S<< "\nD " << D << "\nP " << P << "\nR " <<R << "\nmu " << mu << "\nkappa " << kappa<< "\nkx " << kx << "\nkz " <<kz << "\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "Eres " << Eres<< "\nvresz " << vresz << "\nwtau_sq " << wtau_sq << "\nE_kin " <<  Ekin[p][i] << "\nmu_adiabatic" << eql_dstr[p].M_adiabatic.at(i);
//RK step-1//#################################################################################################################################################################################################################################################################################################
				k1=z_rk(eql_dstr[p].ppar.at(i),electron.mass,gama);
				l1=p_par_rk(eql_dstr[p].pper.at(i),eql_dstr[p].eta.at(i),electron.mass,kz,w_h,dwh_ds,gama, wtau_sq);
				m1=p_per_rk(eql_dstr[p].ppar.at(i),eql_dstr[p].pper.at(i),eql_dstr[p].eta.at(i),electron.mass,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
				n1=eta_rk(electron.mass,Constants::w_wave,kz,eql_dstr[p].ppar.at(i),w_h,gama);
				o1=lamda_rk(eql_dstr[p].ppar.at(i),eql_dstr[p].lamda.at(i),electron.mass,gama);
				p1=alpha_rk(eql_dstr[p].pper.at(i),eql_dstr[p].alpha.at(i),eql_dstr[p].eta.at(i), electron.mass, kz, w_h, Constants::w_wave ,dwh_ds,gama, wtau_sq);
				q1=aeq_rk(eql_dstr[p].ppar.at(i),eql_dstr[p].pper.at(i),eql_dstr[p].alpha.at(i), eql_dstr[p].eta.at(i), eql_dstr[p].aeq.at(i), kappa, gama, Bwc); 
				//std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn " << n1<< "\no1 " << o1 << "\np1 " << p1 << "\nq1 " << q1 <<"\n";	
//RK step-1//#################################################################################################################################################################################################################################################################################################
				Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);
				ns_e = electron.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1); ns_O = oxygen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1); ns_He = helium.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);
				wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
				dwh_ds=dwh_dsf(w_h,eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1);
				p_mag = sqrt((eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1)*(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1)+(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1)*(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1));
		        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(electron.mass*electron.mass*Constants::c*Constants::c*Constants::c*Constants::c))/(electron.mass*Constants::c*Constants::c);
				//Assign no interaction values.
				kz=0;
				wtau_sq=0;
				w1=0;
				w2=0;
				R1=0;
				R2=0;
				beta=0;
				kappa=0;
				Bwc=0;
				if(Constants::interaction)
				{
					wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
					stix = stix_parameters(Constants::w_wave, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
					S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
					disp = dispersion(S,P,R,L,D,Constants::theta0,Constants::w_wave);
					mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
					whistlers(p,i,mu,P,D,S,kz,eql_dstr[p].zeta.at(i)+0.5*Constants::h*k1,time_sim[p][i]+0.5*Constants::h*k1, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
					Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
	    			Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
					Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
				}
//RK step-2//#################################################################################################################################################################################################################################################################################################
				k2=z_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,electron.mass,gama);  
				l2=p_par_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1,electron.mass,kz,w_h,dwh_ds,gama,wtau_sq);
				m2=p_per_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1,electron.mass,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
				n2=eta_rk(electron.mass,Constants::w_wave,kz,eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,w_h,gama);
				o2=lamda_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o1,electron.mass,gama);
				p2=alpha_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1, eql_dstr[p].alpha.at(i)+0.5*(Constants::h)*p1,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1, electron.mass,  kz,  w_h, Constants::w_wave ,dwh_ds,gama, wtau_sq);
				q2=aeq_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l1,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m1,eql_dstr[p].alpha.at(i)+0.5*(Constants::h)*p1, eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n1, eql_dstr[p].aeq.at(i)+0.5*(Constants::h)*q1, kappa, gama, Bwc); 
				//std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\n";
//RK step-2//#################################################################################################################################################################################################################################################################################################
				Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);
				ns_e = electron.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2); ns_O = oxygen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2); ns_He = helium.density(eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);
				wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
				dwh_ds=dwh_dsf(w_h,eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2);
				p_mag = sqrt((eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2)*(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2)+(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2)*(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2));
		        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(electron.mass*electron.mass*Constants::c*Constants::c*Constants::c*Constants::c))/(electron.mass*Constants::c*Constants::c);
				//Assign no interaction values.
				kz=0;
				wtau_sq=0;
				w1=0;
				w2=0;
				R1=0;
				R2=0;
				beta=0;
				kappa=0;
				Bwc=0;
				if(Constants::interaction)
				{
					wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
					stix = stix_parameters(Constants::w_wave, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
					S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
					disp = dispersion(S,P,R,L,D,(Constants::theta0),Constants::w_wave);
					mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
					whistlers(p,i,mu,P,D,S,kz,eql_dstr[p].zeta.at(i)+0.5*Constants::h*k2,time_sim[p][i]*0.5*Constants::h, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
					Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
	    			Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
					Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
				}			
//RK step-3//#################################################################################################################################################################################################################################################################################################
				k3=z_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,electron.mass,gama);  
				l3=p_par_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2,electron.mass,kz,w_h,dwh_ds,gama, wtau_sq);
				m3=p_per_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2,electron.mass,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
				n3=eta_rk(electron.mass,Constants::w_wave,kz,eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,w_h,gama);
				o3=lamda_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].lamda.at(i)+0.5*(Constants::h)*o2,electron.mass,gama);
				p3=alpha_rk(eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2, eql_dstr[p].alpha.at(i)+0.5*(Constants::h)*p2,eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2, electron.mass,  kz,  w_h, Constants::w_wave ,dwh_ds,gama, wtau_sq);
				q3=aeq_rk(eql_dstr[p].ppar.at(i)+0.5*(Constants::h)*l2,eql_dstr[p].pper.at(i)+0.5*(Constants::h)*m2,eql_dstr[p].alpha.at(i)+0.5*(Constants::h)*p2, eql_dstr[p].eta.at(i)+0.5*(Constants::h)*n2, eql_dstr[p].aeq.at(i)+0.5*(Constants::h)*q2, kappa, gama, Bwc); 
				//std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 << "\n";
//RK step-3//#################################################################################################################################################################################################################################################################################################
				Bmag=Bmag_dipole(eql_dstr[p].lamda.at(i)+(Constants::h)*o3);
				ns_e = electron.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3); ns_O = oxygen.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3);  ns_H = hydrogen.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3); ns_He = helium.density(eql_dstr[p].lamda.at(i)+(Constants::h)*o3);
				wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag);
				dwh_ds=dwh_dsf(w_h,eql_dstr[p].lamda.at(i)+(Constants::h)*o3);
				p_mag = sqrt((eql_dstr[p].ppar.at(i)+(Constants::h)*l3)*(eql_dstr[p].ppar.at(i)+(Constants::h)*l3)+(eql_dstr[p].pper.at(i)+(Constants::h)*m3)*(eql_dstr[p].pper.at(i)+(Constants::h)*m3));
		        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(electron.mass*electron.mass*Constants::c*Constants::c*Constants::c*Constants::c))/(electron.mass*Constants::c*Constants::c);
				//Assign no interaction values.
				kz=0;
				wtau_sq=0;
				w1=0;
				w2=0;
				R1=0;
				R2=0;
				beta=0;
				kappa=0;
				Bwc=0;
				if(Constants::interaction)
				{
					wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
					stix = stix_parameters(Constants::w_wave, wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
					S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
					disp = dispersion(S,P,R,L,D,(Constants::theta0),Constants::w_wave);
					mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
					whistlers(p,i,mu,P,D,S,kz,eql_dstr[p].zeta.at(i)+Constants::h*k3,time_sim[p][i]+Constants::h, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
					Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
	    			Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
					Bell_params(electron.mass,electron.charge,eql_dstr[p].ppar.at(i)+(Constants::h)*l3,eql_dstr[p].pper.at(i)+(Constants::h)*m3,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
				}
//RK step-4//#################################################################################################################################################################################################################################################################################################																								
				k4=z_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,electron.mass,gama);  
				l4=p_par_rk(eql_dstr[p].pper.at(i)+(Constants::h)*m3,eql_dstr[p].eta.at(i)+(Constants::h)*n3,electron.mass,kz,w_h,dwh_ds,gama,wtau_sq);
				m4=p_per_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,eql_dstr[p].pper.at(i)+(Constants::h)*m3,eql_dstr[p].eta.at(i)+(Constants::h)*n3,electron.mass,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
				n4=eta_rk(electron.mass,Constants::w_wave,kz,eql_dstr[p].ppar.at(i)+(Constants::h)*l3,w_h,gama);
				o4=lamda_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,eql_dstr[p].lamda.at(i)+(Constants::h)*o3,electron.mass,gama);
				p4=alpha_rk(eql_dstr[p].pper.at(i)+(Constants::h)*m3, eql_dstr[p].alpha.at(i)+(Constants::h)*p3,eql_dstr[p].eta.at(i)+(Constants::h)*n3, electron.mass,  kz,  w_h, Constants::w_wave ,dwh_ds,gama,wtau_sq);
				q4=aeq_rk(eql_dstr[p].ppar.at(i)+(Constants::h)*l3,eql_dstr[p].pper.at(i)+(Constants::h)*m3,eql_dstr[p].alpha.at(i)+(Constants::h)*p3, eql_dstr[p].eta.at(i)+(Constants::h)*n3, eql_dstr[p].aeq.at(i)+(Constants::h)*q3, kappa, gama, Bwc); 
				//std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\n";
//RK step-4//#################################################################################################################################################################################################################################################################################################																								


//RK approx//#################################################################################################################################################################################################################################################################################################
				new_zeta  = eql_dstr[p].zeta.at(i)+((Constants::h)/6)*(k1+2*k2+2*k3+k4); 
				new_ppar  = eql_dstr[p].ppar.at(i)+((Constants::h)/6)*(l1+2*l2+2*l3+l4);
				new_pper  = eql_dstr[p].pper.at(i)+((Constants::h)/6)*(m1+2*m2+2*m3+m4);
				new_eta   = eql_dstr[p].eta.at(i)+((Constants::h)/6)*(n1+2*n2+2*n3+n4);
				new_lamda = eql_dstr[p].lamda.at(i)+((Constants::h)/6)*(o1+2*o2+2*o3+o4);
				new_alpha = eql_dstr[p].alpha.at(i)+((Constants::h)/6)*(p1+2*p2+2*p3+p4);
				new_aeq   = eql_dstr[p].aeq.at(i)+((Constants::h)/6)*(q1+2*q2+2*q3+q4);
				deta_dt[p][i+1]=((Constants::h)/6)*(n1+2*n2+2*n3+n4);
//RK approx//#################################################################################################################################################################################################################################################################################################
		
				//To find minimum-maximum movements in latitude for the particles.
				if(abs(new_lamda - eql_dstr[p].lamda.at(i) ) < min_diff)
				{
					min_diff = abs(new_lamda - eql_dstr[p].lamda.at(i));
				}
				if(abs(new_lamda - eql_dstr[p].lamda.at(i)) > max_diff)
				{
					max_diff = abs(new_lamda - eql_dstr[p].lamda.at(i));
				}

				Beq=Bmag_dipole(0);								//Why call it everytime? Isn't it Beq0 always?
				B_lam=Bmag_dipole(new_lamda);
				
				//salphaeq=sin(new_alpha)*sqrt(Beq/B_lam);		//Not reliable way. It's better to make a runge kutta estimation.
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
				//std::cout<<"\n"<< new_alpha << " " << new_zeta << " " << new_ppar<< " " << new_pper<< " " << new_eta << " " <<new_lamda*Constants::R2D<< " " <<new_aeq ;
				//std::cout<<"\nParticle "<<p<<" lamda: " << new_lamda*Constants::R2D<<" p.a: "<<new_alpha*Constants::R2D<< " upar: " << new_upar << " uper: "<<new_uper;
				
//if particle crosses satellite//#################################################################################################################################################################################################################################################################################################
		        if( ODPT.crossing(new_lamda*Constants::R2D, eql_dstr[p].lamda.at(i-1)*Constants::R2D, Constants::L_shell) )	 
		        {										//i was increased!
					std::cout<<"\nParticle "<< p <<" at: "<<new_lamda*Constants::R2D<< " just crossed the satellite, at: "<< new_time << " simulation seconds\n";
		        	//Store it's state(it's before crossing the satellite!).
		        	ODPT.store(p, eql_dstr[p].lamda.at(i-1), eql_dstr[p].lamda.at(i-1), eql_dstr[p].uper.at(i-1) , eql_dstr[p].upar.at(i-1), eql_dstr[p].ppar.at(i-1), eql_dstr[p].pper.at(i-1), eql_dstr[p].alpha.at(i-1), eql_dstr[p].alpha2.at(i-1), eql_dstr[p].aeq.at(i-1), eql_dstr[p].aeq2.at(i-1), eql_dstr[p].eta.at(i-1), eql_dstr[p].M_adiabatic.at(i-1), eql_dstr[p].time.at(i-1));  	
		        	
		        }
//if particle crosses satellite//#################################################################################################################################################################################################################################################################################################

				//if(eql_dstr[p].lamda.at(i)>0) {	//Stop at equator, 
				//	break;}															
			
			}				
		}

		//std::cout<<"\nMin and max lat move "<<min_diff <<", "<<max_diff<<" respectively\n";
		//std::cout<<"\n"<<ODPT.lamda.at(0)<<" " << ODPT.lamda.at(1);
		auto rk_stop = std::chrono::high_resolution_clock::now();  
		auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(rk_stop - rk_start);
		real rk_time = duration2.count() * pow(10,-6);
		std::cout << "\nExecution time " << std::setprecision(8) << rk_time <<" seconds\n" ;

//--------------------------------------------------------------- RUNGE KUTTA: END -------------------------------------------------------------------//

//--------------------------------------------------------------- CATCH EXCEPTIONS -------------------------------------------------------------------//

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
		else if(exception == 97)
		{
			std::cout<< "\n" << "Caught exception: Not enough RAM.\n";
			return EXIT_FAILURE;
		}
	}
	

	catch(std::bad_alloc& ex) //doesn't work.
	{
		std::cout<< "\n Out of Memory Error.\n";
		return EXIT_FAILURE;
	}

//---------------------------------------------------------- CATCH EXCEPTIONS : END -----------------------------------------------------------------//	


//------------------------------------------------------------ OUTPUT DATA HDF5 ---------------------------------------------------------------------//
	
//-----------------------------------------------FOR ALL PARTICLES---------------------------------------------------//	
	std::vector<std::vector<real>> aeq2_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>>  aeq_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>>  eta_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>>lamda_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>> time_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>> upar_plot(Constants::population, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	//ASSIGN FROM STRUCT INTO 2D VECTORS. BETTER WAY?
	for(int p=0; p<Constants::population; p++)
	{
	    for(unsigned long int i=0; i<eql_dstr[p].lamda.size(); i++)  //Fill 2D vector from vector of structs.
	    {
	        aeq_plot[p][i]  = eql_dstr[p].aeq.at(i);  
	        eta_plot[p][i]  = eql_dstr[p].eta.at(i);  
	        aeq2_plot[p][i] = eql_dstr[p].aeq2.at(i);  
	        lamda_plot[p][i]= eql_dstr[p].lamda.at(i);  
	    	time_plot[p][i] = eql_dstr[p].time.at(i);
	    	upar_plot[p][i] = eql_dstr[p].upar.at(i);
	    }
	}
	//Create hdf5 file.
	h5::File file1("particles.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	//Create datasets and write vectors to them(shortcut syntax).
	h5::DataSet dataset1_lamda =  file1.createDataSet("lamda 2D",       	lamda_plot);
	h5::DataSet dataset1_aeq   =  file1.createDataSet("aeq 2D",         	aeq_plot);
	h5::DataSet dataset1_aeq2  =  file1.createDataSet("aeq2 2D",        	aeq2_plot);
	h5::DataSet dataset1_eta   =  file1.createDataSet("eta 2D",         	eta_plot);
	h5::DataSet dataset1_time  =  file1.createDataSet("simulation time",	time_plot);
	h5::DataSet dataset1_upar  =  file1.createDataSet("upar",	    		upar_plot);

//--------------------------------------------FOR DETECTED PARTICLES------------------------------------------------//	
	//Create hdf5 file for particles that cross the satellite.
	h5::File file2("detected.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);

	h5::DataSet dataset2_lamda = file2.createDataSet("detected_lamda 1D", ODPT.lamda);
	h5::DataSet dataset2_time  = file2.createDataSet("detected_time 1D",  ODPT.time);
	h5::DataSet dataset2_id    = file2.createDataSet("detected_id 1D",  ODPT.id);
	h5::DataSet dataset2_aeq   = file2.createDataSet("detected_aeq 1D", ODPT.aeq);
	h5::DataSet dataset2_alpha = file2.createDataSet("detected_alpha 1D", ODPT.alpha);
	h5::DataSet telescope_lamda= file2.createDataSet("telescope_lamda scalar", ODPT.latitude);

//----------------------------------------------------------- OUTPUT DATA HDF5 : END -------------------------------------------------------------//

return 0; 

}


