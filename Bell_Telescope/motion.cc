#include "headers/motion.h"

void motion(int p, real lamda, real alpha, real aeq, real ppar, real pper, real upar, real uper, real zeta, real M_adiabatic, real eta, real time, Telescope &ODPT)
{
    //Declare function's variables. Once for each particle. When parallel, declare xcore times?
    real new_lamdae;
    real ns_e, wc_e, wps_e, ns_O, wc_O, wps_O ,ns_H, wc_H, wps_H, ns_He, wc_He, wps_He, w_h;
    real Bmag;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;
    real gama,w1,w2,R1,R2,beta,wtau_sq;
    real S,D,P,R,L,mu,dwh_ds,vresz,Eres,kappa,kx,kz;
    real Bxwc,Bzwc,Bywc,Exwc,Eywc,Ezwc,Ewc,Bwc;
    real p_mag;
    int i=0;
    //Tuples
    std::tuple<real, real, real, real, real> stix;
    std::tuple<real, real, real, real> disp;

    //Objects for each specie. Used inside "motion" function. //Change that, it's called once for each particle...
    Species electron(Constants::m_e,  Constants::q_e, 1); 
    Species oxygen  (Constants::m_O,  Constants::q_i, 0.006); 
    Species hydrogen(Constants::m_H,  Constants::q_i, 0.94); 
    Species helium  (Constants::m_He, Constants::q_i, 0.054);


    std::vector <std::vector<real>> time_sim  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //Optional

    //Make vectors to save data throughout the Runge Kutta only if needed.(0) 
    //Create vectors instead of big arrays to allocate in heap.
    //Each row(population) has a vector which has cols(Nsteps+1) number of elements, each element initialized to 0.

    //std::vector <std::vector<real>> Ekin      (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) ); 
    //std::vector <std::vector<real>> Exw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> Eyw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> Ezw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> Bxw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> Byw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> Bzw_out   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> Bw_out    (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> Ew_out    (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> kappa_out (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> kx_out	  (Constants::population, std=::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> kz_out	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> L_stix	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> S_stix	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> D_stix	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> P_stix	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> R_stix	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> wh_out    (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> mu_out	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> vresz_o   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> Eres_o	  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> gama_out  (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> deta_dt   (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>>dwh_dt_out (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>>B_earth_out(Constants::population, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>>Phi_out    (Constants::population, std::vector<real> (Constants::Nsteps + 1, 0);
    
	//std::cout.precision(32);			//Output 16 decimal precise
	//std::cout<<std::scientific;		//For e notation representation
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
            whistlers(p,i,mu,P,D,S,kz,zeta,time_sim[p][i], Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
            //Total fields.
            Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
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
        //std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He << "\nlamda " << lamda<< "\naeq " << eql_dstr[p].aeq2.at(i) << "\nBxw " << Bxwc << "\nByw "<<Bywc<< "\nBzw "<<Bzwc<< "\nBw " << Bw_out[p][i] << "\nExw " << Exwc<< "\nEyw " << Eywc << "\nEzw " << Ezwc << "\nEw " << Ew_out[p][i] << "\nL " << L << "\nS " <<S<< "\nD " << D << "\nP " << P << "\nR " <<R << "\nmu " << mu << "\nkappa " << kappa<< "\nkx " << kx << "\nkz " <<kz << "\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "Eres " << Eres<< "\nvresz " << vresz << "\nwtau_sq " << wtau_sq << "\nE_kin " <<  Ekin[p][i] << "\nmu_adiabatic" << eql_dstr[p].M_adiabatic.at(i);
        
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
            stix = stix_parameters(wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
            S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
            disp = dispersion(S,P,R,L,D);
            mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
            whistlers(p,i,mu,P,D,S,kz,zeta+0.5*Constants::h*k1,time_sim[p][i]+0.5*Constants::h*k1, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
            Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
            Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
            Bell_params(ppar+0.5*(Constants::h)*l1,pper+0.5*(Constants::h)*m1,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
        }
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
            stix = stix_parameters(wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
            S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
            disp = dispersion(S,P,R,L,D);
            mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
            whistlers(p,i,mu,P,D,S,kz,zeta+0.5*Constants::h*k2,time_sim[p][i]*0.5*Constants::h, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
            Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
            Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
            Bell_params(ppar+0.5*(Constants::h)*l2,pper+0.5*(Constants::h)*m2,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
        }			
        
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
            stix = stix_parameters(wc_e, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
            S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
            disp = dispersion(S,P,R,L,D);
            mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
            whistlers(p,i,mu,P,D,S,kz,zeta+Constants::h*k3,time_sim[p][i]+Constants::h, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
            Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
            Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
            Bell_params(ppar+(Constants::h)*l3,pper+(Constants::h)*m3,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
        }
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
            std::cout<<"\n\nParticle "<<p<<" with aeq0="<<eql_dstr[p].aeq[0]*Constants::R2D<<" breaks.";
            //std::cout<<"\n"<< alpha << " " << zeta << " " << ppar<< " " << pper<< " " << eta << " " <<lamda<< " " <<aeq ;
            //std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He << "\nlamda " <<"\nBwc " << Bwc << "\nEwc "<< Ewc << "\nL " << L << "\nS " <<S<< "\nD " << D << "\nP " << P << "\nR " <<R << "\nmu " << mu << "\nkappa " << kappa<< "\nkx " << kx << "\nkz " <<kz << "\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "Eres " << Eres<< "\nvresz " << vresz << "\nwtau_sq " << wtau_sq << "\nmu_adiabatic" << M_adiabatic;
                    
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
        upar  =  ppar   /  (Constants::m_e*gama);
        uper  =  pper   /  (Constants::m_e*gama);

        deta_dt[p][i+1] = (Constants::h/6)*(n1+2*n2+2*n3+n4);

        //B_lam    =  Bmag_dipole(lamda);    
        //M_adiabatic = (pper*pper)/(2*Constants::m_e*B_lam); 

        //Go to next timestep
        time  = time + Constants::h; 
        
        eql_dstr[p].update_state(aeq, time);

        //std::cout<<"\n\nzeta "<< zeta << "\nppar "<< ppar<< "\npper " << pper<< "\neta " << eta << "\nlamda " <<lamda<< "\nalpha "<< alpha << "\naeq " <<aeq ;

        i++;  

        //Stop at equator
        //if(eql_dstr[p].lamda.at(i)>0) {	
        //	break;}	
    }
}