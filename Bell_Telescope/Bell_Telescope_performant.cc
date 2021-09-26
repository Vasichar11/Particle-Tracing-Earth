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
	Species electron(Constants::m_e,  Constants::q_e, 1); 
	Species oxygen  (Constants::m_O,  Constants::q_i, 0.006); 
	Species hydrogen(Constants::m_H,  Constants::q_i, 0.94); 
	Species helium  (Constants::m_He, Constants::q_i, 0.054);	
 
	//Position of the Particle Telescope.		
	Telescope ODPT(Constants::telescope_lamda, Constants::L_shell);		

//-------------------------------------------------------------DISTRIBUTION OF PARTICLES----------------------------------------------------------------//
	//Object for particles.
	std::cout<<"\n\nParticle testing population: " << Constants::test_pop << "\n\nWave interaction: "<<Constants::interaction*Constants::By_wave<< " T" ;
	std::cout<<"\n\nEta distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::eta_start_d << "  "<< " " << Constants::eta_end_d <<"|\n";
	std::cout<<"\nWith aeq distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::aeq_start_d << "  "<< " " << Constants::aeq_end_d <<"|\n";
	std::cout<<"\nWith lamda distribution in degrees"<<"\n|From "<<" To|";
	std::cout<<"\n| "<<Constants::lamda_start_d << "  "<< " " << Constants::lamda_end_d <<"|\n";

	Particles single; //Single particle struct.
	std::vector<Particles> eql_dstr(Constants::test_pop, single);	//Vector of structs for particle distribution.
						//equally distributed
	
	real eta0,aeq0,lamda0;
	real Blam0,salpha0,alpha0;
	real Beq0 = Bmag_dipole(0);	//Beq isn't always Beq0?

	for(int e=0, p=0; e<Constants::eta_dstr; e++)		
	{   																							 
		eta0 = (Constants::eta_start_d + e*Constants::eta_step_d) * Constants::D2R;			   														
		
		for(int a=0; a<Constants::aeq_dstr; a++)											
		{	
			aeq0 = (Constants::aeq_start_d + a*Constants::aeq_step_d) * Constants::D2R;	         	  	    

			for(int l=0; l<Constants::lamda_dstr; l++, p++)
		    {
		    	lamda0 = (Constants::lamda_start_d + l*Constants::lamda_step_d) * Constants::D2R;  	
				
				//Find P.A at lamda0.
				Blam0=Bmag_dipole(lamda0);
				salpha0=sin(aeq0)*sqrt(Blam0/Beq0); //(2.20) Bortnik thesis
				if( (salpha0<-1) || (salpha0>1) || (salpha0==0) ) {eql_dstr.pop_back(); p--; continue; } //Exclude these particles.
				alpha0=asin(salpha0);		//If aeq0=150 => alpha0=arcsin(sin(150))=30 for particle in equator.Distribute in alpha instead of aeq?		 	
				eql_dstr[p].initialize(eta0,aeq0,alpha0,lamda0,Constants::Ekev0,Blam0,0,0);

				//Print initial state of particles.
				//std::cout<<"\nParticle"<<p<<" with eta0: "<< eql_dstr[p].eta.at(0)*Constants::R2D <<", aeq0: "<< aeq0*Constants::R2D <<" gives alpha0: "<<alpha0*Constants::R2D<<" and lamda0: "<<lamda0*Constants::R2D<<"\n";				
			}
		}	
	}
	int64_t track_pop = eql_dstr.size(); //Population of particles that will be tracked.
	std::cout<<"\n"<<Constants::test_pop - track_pop<<" Particles were excluded from the initial population due to domain issues.\nThe particle population for the tracer is now: "<<track_pop<<"\n";
//-------------------------------------------------------------DISTRIBUTION OF PARTICLES:END------------------------------------------------------------//

	std::cout.precision(8);			//Output 16 decimal precise
	std::cout<<std::scientific;		//For e notation representation

//--------------------------------------------------------------------RUNGE KUTTA-----------------------------------------------------------------------//

	//Temp variables for Runge Kutta
	real lamda , zeta, uper , upar, ppar, pper, alpha, aeq, eta, time;
	real new_lamda;

	real ns_e, wc_e, wps_e, ns_O, wc_O, wps_O ,ns_H, wc_H, wps_H, ns_He, wc_He, wps_He, w_h;
	real Bmag;
	real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;
	real gama,w1,w2,R1,R2,beta,wtau_sq;
	real S,D,P,R,L,mu,dwh_ds,vresz,Eres,kappa,kx,kz;
	real Bxwc,Bzwc,Bywc,Exwc,Eywc,Ezwc,Bwc;
	real p_mag;
	std::tuple<real, real, real, real, real> stix;
	std::tuple<real, real, real, real> disp;

	
	std::vector <std::vector<real>> time_sim  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );


	//Temp check size of vectors and if it can be allocated. Overheads are not included.
	int64_t szv = ( (time_sim.size() * time_sim[0].size() * sizeof(real)) + (time_sim.size()*sizeof(time_sim[0])) + sizeof(time_sim) ) * 40 ; 
	std::cout<< "\nTo allocate in memory all simulation's 2D Vectors you will need more than: ~ " << szv <<"B * 40" << szv*pow(10,-9) << " free GB \n";
	
	//bounce is ~1s, therefore detection every 0.1s => simulation time/0.1 instead of vectors of Nsteps size.
	int64_t szv_performant = ((time_sim.size() * (Constants::t/1) * sizeof(real)) + (time_sim.size()*sizeof(time_sim[0])) + sizeof(time_sim) ) * 10 ;
	
	std::cout<< "\nTo allocate in memory only detected particles you will need more than: ~ " << szv_performant <<"B * 10 = " << szv_performant*pow(10,-9) << " free GB \n";


	//Make vectors to save data throughout the Runge Kutta only if needed.(0) 

	//Initializations of arrays(Attributes of Particles are vectors inside Particles struct).
	//Create vectors instead of big arrays to allocate in heap.
	//Each row(population) has a vector which has cols(Nsteps+1) number of elements, each element initialized to 0.
	
	//std::vector <std::vector<real>> Ekin      (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) ); 
	//std::vector <std::vector<real>> Exw_out   (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> Eyw_out   (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> Ezw_out   (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> Bxw_out   (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> Byw_out   (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> Bzw_out   (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> Bw_out    (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> Ew_out    (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> kappa_out (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> kx_out	  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> kz_out	  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> L_stix	  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> S_stix	  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> D_stix	  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> P_stix	  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> R_stix	  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> wh_out    (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> mu_out	  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> vresz_o   (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> Eres_o	  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> gama_out  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>> deta_dt   (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>>dwh_dt_out (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>>B_earth_out(track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
	//std::vector <std::vector<real>>Phi_out    (track_pop, std::vector<real> (Constants::Nsteps + 1, 0);


	std::cout<<"\n";
	
	auto rk_start = std::chrono::high_resolution_clock::now();

	int i;
	int broken=0;
	
	//Assign no interaction values.
    kz=wtau_sq=w1=w2=R1=R2=beta=kappa=Bwc=0;
	
	for(int p=0; p<track_pop; p++)     //Loop for all particles
	{
		
		i=0;
		lamda = eql_dstr[p].lamda.at(i);
		alpha = eql_dstr[p].alpha.at(i);
		aeq   = eql_dstr[p].aeq.at(i);
		ppar  = eql_dstr[p].ppar.at(i);
		pper  = eql_dstr[p].pper.at(i);
		upar  = eql_dstr[p].upar.at(i);
		uper  = eql_dstr[p].uper.at(i);
		eta   = eql_dstr[p].eta.at(i);
		zeta   = eql_dstr[p].zeta.at(i);
		time  = eql_dstr[p].time.at(i);
		//M_adiabatic = eql_dstr[p].M_adiabatic.at(i);
		
		
		std::cout<<"\rParticle "<<p<<" is bouncing"<<std::flush;
		
		while(i<Constants::Nsteps) 
		{
			
			Bmag=Bmag_dipole(lamda);
			//Densities.
			ns_e = electron.density(lamda); ns_O = oxygen.density(lamda);  ns_H = hydrogen.density(lamda); ns_He = helium.density(lamda);
			//Cyclotron frequencies.
			wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);  wc_He = helium.wc(Bmag);
			//dwh_dt.
			dwh_ds=dwh_dsf(w_h,lamda);
			//From bell parameters function(temp for no interaction).
			p_mag = sqrt(ppar*ppar+pper*pper);
			gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
			


			//Keep these values in arrays only if needed.(1)
			//B_earth_out[p][i]=Bmag; 
			//wh_out[p][i] = w_h;
			//dwh_dt_out[p][i] = dwh_ds;
			//gama_out[p][i]=gama;
		

			if(Constants::interaction)
			{
				//Plasma frequencies squared.
				wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
				//Stix.
				stix = stix_parameters(wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
				S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix); //<get> efficiency comperable to accessing a member of a struct
				//Dispersion.
				disp = dispersion(S,P,R,L,D);
				mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
				//Fields.
				whistlers(p,i,mu,P,D,S,kz,zeta,time, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
				//Total fields.
				//Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
				Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);

				//Bell_params.
				Bell_params(ppar,pper,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
				//vresonant.
				vres_f(kz,w_h,alpha,vresz,Eres); //vres?
				
				//Keep these values in arrays only if needed.(2)
				//B_earth_out[p][i]=Bmag; 
				//vresz_o[p][i]= vresz;
				//Eres_o[p][i] = Eres;
				//Exw_out[p][i]= Exwc;
				//Eyw_out[p][i]= Eywc;
				//Ezw_out[p][i]= Ezwc;
				//Bxw_out[p][i]= Bxwc;
				//Byw_out[p][i]= Bywc;
				//Bzw_out[p][i]= Bzwc;
				//S_stix[p][i] = S;
				//D_stix[p][i] = D;
				//P_stix[p][i] = P;
				//R_stix[p][i] = R;
				//L_stix[p][i] = L;
				//Ew_out[p][i] = Ewc;
				//Bw_out[p][i] = Bwc;
				//kx_out[p][i] = kx;
				//kz_out[p][i] = kz;
				//mu_out[p][i] = mu;
				//kappa_out[p][i]  = kappa;
				//wh_out[p][i] = w_h;
				//dwh_dt_out[p][i] = dwh_ds;
				//gama_out[p][i]   = gama;
				//Phi_out[p][i]    = Ph_w;
				//Ekin[p][i]=((gama-1)*Constants::m_e*Constants::c*Constants::c)*6.2415e15; //Joules in kev
			
			}
			//Check print parameters
			//std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He << "\nBwc " << Bwc <<"\nBxwc "<<Bxwc <<"\nBywc "<<Bywc<<"\nBzwc "<<Bzwc<< "\nEwc " << Ewc<<"\nExwc "<<Exwc<<"\nEywc "<<Eywc<<"\nEzwc "<<Ezwc<<"\nL " << L << "\nS " <<S<< "\nD " << D << "\nP " << P << "\nR " <<R << "\nmu " << mu << "\nkappa " << kappa<< "\nkx " << kx << "\nkz " <<kz << "\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "Eres " << Eres<< "\nvresz " << vresz << "\nwtau_sq " << wtau_sq << "\nmu_adiabatic" << M_adiabatic;
			
			//RK step-1//#################################################################################################################################################################################################################################################################################################
			k1=z_rk(ppar,gama);
			l1=p_par_rk(pper,eta,kz,w_h,dwh_ds,gama, wtau_sq);
			m1=p_per_rk(ppar,pper,eta,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
			n1=eta_rk(kz,ppar,w_h,gama);
			o1=lamda_rk(ppar,lamda,gama);
			p1=alpha_rk(pper,alpha,eta, kz, w_h,dwh_ds,gama, wtau_sq);
			q1=aeq_rk(ppar,pper,alpha, eta, aeq, kappa, gama, Bwc); 
			//std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn " << n1<< "\no1 " << o1 << "\np1 " << p1 << "\nq1 " << q1 <<"\n";	
			
			Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o1);
			ns_e = electron.density(lamda+0.5*(Constants::h)*o1); ns_O = oxygen.density(lamda+0.5*(Constants::h)*o1);  ns_H = hydrogen.density(lamda+0.5*(Constants::h)*o1); ns_He = helium.density(lamda+0.5*(Constants::h)*o1);
			wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
			dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o1);
			p_mag = sqrt((ppar+0.5*(Constants::h)*l1)*(ppar+0.5*(Constants::h)*l1)+(pper+0.5*(Constants::h)*m1)*(pper+0.5*(Constants::h)*m1));
			gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
			

			if(Constants::interaction)
			{
				wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
				stix = stix_parameters(wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
				S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
				disp = dispersion(S,P,R,L,D);
				mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
				whistlers(p,i,mu,P,D,S,kz,zeta+0.5*Constants::h*k1,time+0.5*Constants::h*k1, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
				//Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
				Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
				Bell_params(ppar+0.5*(Constants::h)*l1,pper+0.5*(Constants::h)*m1,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
			}
			//std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He << "\nBwc " << Bwc <<"\nBxwc "<<Bxwc <<"\nBywc "<<Bywc<<"\nBzwc "<<Bzwc<< "\nEwc " << Ewc<<"\nExwc "<<Exwc<<"\nEywc "<<Eywc<<"\nEzwc "<<Ezwc<<"\nL " << L << "\nS " <<S<< "\nD " << D << "\nP " << P << "\nR " <<R << "\nmu " << mu << "\nkappa " << kappa<< "\nkx " << kx << "\nkz " <<kz << "\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "Eres " << Eres<< "\nvresz " << vresz << "\nwtau_sq " << wtau_sq << "\nmu_adiabatic" << M_adiabatic;
			
			//RK step-2//#################################################################################################################################################################################################################################################################################################
			k2=z_rk(ppar+0.5*(Constants::h)*l1,gama);  
			l2=p_par_rk(pper+0.5*(Constants::h)*m1,eta+0.5*(Constants::h)*n1,kz,w_h,dwh_ds,gama,wtau_sq);
			m2=p_per_rk(ppar+0.5*(Constants::h)*l1,pper+0.5*(Constants::h)*m1,eta+0.5*(Constants::h)*n1,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
			n2=eta_rk(kz,ppar+0.5*(Constants::h)*l1,w_h,gama);
			o2=lamda_rk(ppar+0.5*(Constants::h)*l1,lamda+0.5*(Constants::h)*o1,gama);
			p2=alpha_rk(pper+0.5*(Constants::h)*m1, alpha+0.5*(Constants::h)*p1,eta+0.5*(Constants::h)*n1,   kz,  w_h, dwh_ds,gama, wtau_sq);
			q2=aeq_rk(ppar+0.5*(Constants::h)*l1,pper+0.5*(Constants::h)*m1,alpha+0.5*(Constants::h)*p1, eta+0.5*(Constants::h)*n1, aeq+0.5*(Constants::h)*q1, kappa, gama, Bwc); 
			//std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\n";
			
			Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o2);
			ns_e = electron.density(lamda+0.5*(Constants::h)*o2); ns_O = oxygen.density(lamda+0.5*(Constants::h)*o2);  ns_H = hydrogen.density(lamda+0.5*(Constants::h)*o2); ns_He = helium.density(lamda+0.5*(Constants::h)*o2);
			wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
			dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o2);
			p_mag = sqrt((ppar+0.5*(Constants::h)*l2)*(ppar+0.5*(Constants::h)*l2)+(pper+0.5*(Constants::h)*m2)*(pper+0.5*(Constants::h)*m2));
			gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
			

			if(Constants::interaction)
			{
				wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
				stix = stix_parameters(wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
				S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
				disp = dispersion(S,P,R,L,D);
				mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
				whistlers(p,i,mu,P,D,S,kz,zeta+0.5*Constants::h*k2,time*0.5*Constants::h, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
				//Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
				Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
				Bell_params(ppar+0.5*(Constants::h)*l2,pper+0.5*(Constants::h)*m2,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
			}			
			//std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He << "\nBwc " << Bwc <<"\nBxwc "<<Bxwc <<"\nBywc "<<Bywc<<"\nBzwc "<<Bzwc<< "\nEwc " << Ewc<<"\nExwc "<<Exwc<<"\nEywc "<<Eywc<<"\nEzwc "<<Ezwc<<"\nL " << L << "\nS " <<S<< "\nD " << D << "\nP " << P << "\nR " <<R << "\nmu " << mu << "\nkappa " << kappa<< "\nkx " << kx << "\nkz " <<kz << "\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "Eres " << Eres<< "\nvresz " << vresz << "\nwtau_sq " << wtau_sq << "\nmu_adiabatic" << M_adiabatic;
			
			//RK step-3//#################################################################################################################################################################################################################################################################################################
			k3=z_rk(ppar+0.5*(Constants::h)*l2,gama);  
			l3=p_par_rk(pper+0.5*(Constants::h)*m2,eta+0.5*(Constants::h)*n2,kz,w_h,dwh_ds,gama, wtau_sq);
			m3=p_per_rk(ppar+0.5*(Constants::h)*l2,pper+0.5*(Constants::h)*m2,eta+0.5*(Constants::h)*n2,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
			n3=eta_rk(kz,ppar+0.5*(Constants::h)*l2,w_h,gama);
			o3=lamda_rk(ppar+0.5*(Constants::h)*l2,lamda+0.5*(Constants::h)*o2,gama);
			p3=alpha_rk(pper+0.5*(Constants::h)*m2, alpha+0.5*(Constants::h)*p2,eta+0.5*(Constants::h)*n2,   kz,  w_h,dwh_ds,gama, wtau_sq);
			q3=aeq_rk(ppar+0.5*(Constants::h)*l2,pper+0.5*(Constants::h)*m2,alpha+0.5*(Constants::h)*p2, eta+0.5*(Constants::h)*n2, aeq+0.5*(Constants::h)*q2, kappa, gama, Bwc); 
			//std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 << "\n";

			Bmag=Bmag_dipole(lamda+(Constants::h)*o3);
			ns_e = electron.density(lamda+(Constants::h)*o3); ns_O = oxygen.density(lamda+(Constants::h)*o3);  ns_H = hydrogen.density(lamda+(Constants::h)*o3); ns_He = helium.density(lamda+(Constants::h)*o3);
			wc_e = w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag);
			dwh_ds=dwh_dsf(w_h,lamda+(Constants::h)*o3);
			p_mag = sqrt((ppar+(Constants::h)*l3)*(ppar+(Constants::h)*l3)+(pper+(Constants::h)*m3)*(pper+(Constants::h)*m3));
			gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
			

			if(Constants::interaction)
			{
				wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
				stix = stix_parameters(wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
				S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
				disp = dispersion(S,P,R,L,D);
				mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
				whistlers(p,i,mu,P,D,S,kz,zeta+Constants::h*k3,time+Constants::h, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
				//Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
				Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
				Bell_params(ppar+(Constants::h)*l3,pper+(Constants::h)*m3,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
			}
			//std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He << "\nBwc " << Bwc <<"\nBxwc "<<Bxwc <<"\nBywc "<<Bywc<<"\nBzwc "<<Bzwc<< "\nEwc " << Ewc<<"\nExwc "<<Exwc<<"\nEywc "<<Eywc<<"\nEzwc "<<Ezwc<<"\nL " << L << "\nS " <<S<< "\nD " << D << "\nP " << P << "\nR " <<R << "\nmu " << mu << "\nkappa " << kappa<< "\nkx " << kx << "\nkz " <<kz << "\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "Eres " << Eres<< "\nvresz " << vresz << "\nwtau_sq " << wtau_sq << "\nmu_adiabatic" << M_adiabatic;
			
			//RK step-4//#################################################################################################################################################################################################################################################################################################																								
			k4=z_rk(ppar+(Constants::h)*l3,gama);  
			l4=p_par_rk(pper+(Constants::h)*m3,eta+(Constants::h)*n3,kz,w_h,dwh_ds,gama,wtau_sq);
			m4=p_per_rk(ppar+(Constants::h)*l3,pper+(Constants::h)*m3,eta+(Constants::h)*n3,w_h,dwh_ds,gama,w1,w2,R1,R2,beta);
			n4=eta_rk(kz,ppar+(Constants::h)*l3,w_h,gama);
			o4=lamda_rk(ppar+(Constants::h)*l3,lamda+(Constants::h)*o3,gama);
			p4=alpha_rk(pper+(Constants::h)*m3, alpha+(Constants::h)*p3,eta+(Constants::h)*n3,   kz,  w_h,dwh_ds,gama,wtau_sq);
			q4=aeq_rk(ppar+(Constants::h)*l3,pper+(Constants::h)*m3,alpha+(Constants::h)*p3, eta+(Constants::h)*n3, aeq+(Constants::h)*q3, kappa, gama, Bwc); 
			//std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\n";



			//Approximate new lamda
			new_lamda = lamda + ((Constants::h)/6)*(o1+2*o2+2*o3+o4);
			
			//Check if NAN to break this particle. Why NAN ?
			if(std::isnan(new_lamda))
			{
				
				std::cout<<"\nParticle "<<p<<" with aeq0="<<eql_dstr[p].aeq[0]*Constants::R2D<<" breaks.\n";
				//std::cout<<"\n"<< alpha << " " << zeta << " " << ppar<< " " << pper<< " " << eta << " " <<lamda<< " " <<aeq ;
				//std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He << "\nlamda " <<"\nBwc " << Bwc << "\nEwc "<< Ewc << "\nL " << L << "\nS " <<S<< "\nD " << D << "\nP " << P << "\nR " <<R << "\nmu " << mu << "\nkappa " << kappa<< "\nkx " << kx << "\nkz " <<kz << "\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "Eres " << Eres<< "\nvresz " << vresz << "\nwtau_sq " << wtau_sq << "\nmu_adiabatic" << M_adiabatic;
				broken++;		
				break;
			}

			//Check crossing. First estimate new latitude. 
			if( ODPT.crossing(new_lamda*Constants::R2D, lamda*Constants::R2D, Constants::L_shell) )	 
			{										
				//std::cout<<"\nParticle "<< p <<" at: "<<new_lamda*Constants::R2D<< " is about to cross the satellite, at: "<< time << " simulation seconds\n";
				//Store its state(it's before crossing the satellite!).
				ODPT.store( p, lamda, uper , upar, alpha, aeq, eta, time);  			        	
			}
			

			//Now approximate all values of Runge Kutta's block.
			lamda =  new_lamda;
			zeta  =  zeta   +  (Constants::h/6)*(k1+2*k2+2*k3+k4);
			ppar  =  ppar   +  (Constants::h/6)*(l1+2*l2+2*l3+l4);
			pper  =  pper   +  (Constants::h/6)*(m1+2*m2+2*m3+m4);
			eta   =  eta    +  (Constants::h/6)*(n1+2*n2+2*n3+n4);
			alpha =  alpha  +  (Constants::h/6)*(p1+2*p2+2*p3+p4);
			aeq   =  aeq    +  (Constants::h/6)*(q1+2*q2+2*q3+q4);
			upar  =  ppar   /  (Constants::m_e*gama); //quite different from zeta?
			uper  =  pper   /  (Constants::m_e*gama);

			//deta_dt[p][i+1] = (Constants::h/6)*(n1+2*n2+2*n3+n4);

			//B_lam    =  Bmag_dipole(lamda);    
			//M_adiabatic = (pper*pper)/(2*Constants::m_e*B_lam); 

			//Go to next timestep
			time  = time + Constants::h; 
			
			//eql_dstr[p].save_state(aeq, alpha, lamda, time);

			//std::cout<<"\n\nzeta "<< zeta << "\nppar "<< ppar<< "\npper " << pper<< "\neta " << eta << "\nlamda " <<lamda<< "\nalpha "<< alpha << "\naeq " <<aeq ;

			i++;  

			//Stop at equator
			//if(eql_dstr[p].lamda.at(i)>0) {	
			//	break;}																

		}				
	}
	//std::cout<<"\nParticles that went out of domain, inside the simulation: "<<broken;
	//std::cout<<"\nMin and max lat move "<<min_diff <<", "<<max_diff<<" respectively\n";
	//std::cout<<"\n"<<ODPT.lamda.at(0)<<" " << ODPT.lamda.at(1);
	auto rk_stop = std::chrono::high_resolution_clock::now();  
	auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(rk_stop - rk_start);
	real rk_time = duration2.count() * pow(10,-6);
	std::cout << "\nExecution time " << std::setprecision(8) << rk_time <<" seconds\n" ;
	
//--------------------------------------------------------------- RUNGE KUTTA: END -------------------------------------------------------------------//




//------------------------------------------------------------ OUTPUT DATA HDF5 ---------------------------------------------------------------------//
/*	
	std::vector<std::vector<real>>  aeq_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>> lamda_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>> time_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	std::vector<std::vector<real>> alpha_plot(track_pop, std::vector<real> (Constants::Nsteps + 1 ,0 ) );
	

	//Assign from struct to 2d vectors.
	for(int p=0; p<track_pop; p++)
	{
	    for(size_t i=0; i<eql_dstr[p].alpha.size(); i++)  
	    {
	    	time_plot[p][i] = eql_dstr[p].time.at(i);
	    	aeq_plot[p][i] = eql_dstr[p].aeq.at(i);
	    	alpha_plot[p][i] = eql_dstr[p].alpha.at(i);
	    	lamda_plot[p][i] = eql_dstr[p].lamda.at(i);

	    }
	}
*/
	h5::File file("h5files/detected.h5", h5::File::ReadWrite | h5::File::Create | h5::File::Truncate);
	
	//Detected particles
	h5::DataSet dataset_lamda      = file.createDataSet("ODPT.lamda", ODPT.lamda);
	h5::DataSet dataset_time       = file.createDataSet("ODPT.time", ODPT.time);
	h5::DataSet dataset_id         = file.createDataSet("ODPT.id", ODPT.id);
	h5::DataSet dataset_alpha      = file.createDataSet("ODPT.alpha", ODPT.alpha);

	//Simulation data and Telescope specification - Scalars 
	h5::DataSet wave_magnitude     = file.createDataSet("By_wave",Constants::By_wave);
	h5::DataSet telescope_lamda    = file.createDataSet("ODPT.latitude", ODPT.latitude);
	h5::DataSet population         = file.createDataSet("population", track_pop);
	h5::DataSet initial_lamda      = file.createDataSet("lamda0", Constants::lamda0);
	h5::DataSet simulation_time    = file.createDataSet("t", Constants::t);
	
	//Saved Particles
	//h5::DataSet dataset_lamda_saved  = file.createDataSet("lamda_plot", lamda_plot);
	//h5::DataSet dataset_alpha_saved  = file.createDataSet("alpha_plot", alpha_plot);
	//h5::DataSet dataset_deta_saved   = file.createDataSet("deta_dt", deta_dt);
	//h5::DataSet dataset_aeq_saved    = file.createDataSet("aeq_plot", aeq_plot);
	//h5::DataSet dataset_time_saved   = file.createDataSet("time_plot", time_plot);




//----------------------------------------------------------- OUTPUT DATA HDF5 : END -------------------------------------------------------------//

return 0; 

}


