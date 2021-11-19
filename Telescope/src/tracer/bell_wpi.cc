#include "headers/bell_wpi.h"

//For wave-particle interaction
void wpi(int p, Particles &single, Telescope &ODPT)
{
	

    real lamda    =  single.lamda.at(0);
    real zeta     =  single.zeta.at(0); 
    real ppar     =  single.ppar.at(0); 
    real pper     =  single.pper.at(0); 
    real eta      =  single.eta.at(0); 
    real alpha    =  single.alpha.at(0); 
    real aeq      =  single.aeq.at(0); 
    //real upar     =  single.upar.at(0); 
    //real uper     =  single.uper.at(0);
    //real deta_dt  =  single.deta_dt.at(0);
    //real Ekin     =  single.Ekin.at(0);
    //real M_adiabatic = single.M_adiabatic.at(0);
    real time     =  single.time.at(0);

    //Declare function's variables. Once for each particle. When parallel, declare xcore times?
    real new_lamda;
    real ns_e,w_h, wps_e, ns_O, wc_O, wps_O ,ns_H, wc_H, wps_H, ns_He, wc_He, wps_He;
    real Bmag;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;
    real gama,w1,w2,R1,R2,beta,wtau_sq;
    real S,D,P,R,L,mu,dwh_ds,vresz,Eres,kappa,kx,kz;
    real Bxwc,Bzwc,Bywc,Exwc,Eywc,Ezwc,Bwc;
    real p_mag;
    //Tuples
    std::tuple<real, real, real, real, real> stix;
    std::tuple<real, real, real, real> disp;
    //Objects for each specie. Used inside "motion" function. //Change that, it's called once for each particle...
    Species electron(Constants::m_e,  Constants::q_e, 1); 
    Species oxygen  (Constants::m_O,  Constants::q_i, 0.006); 
    Species hydrogen(Constants::m_H,  Constants::q_i, 0.94); 
    Species helium  (Constants::m_He, Constants::q_i, 0.054);
    

    
	std::cout.precision(8);			//Output 16 decimal precise
	std::cout<<std::scientific;		//For e notation representation
    
    int i=0;

    while(i<Constants::Nsteps) 
    {
        Bmag=Bmag_dipole(lamda);
        ns_e = electron.density(lamda); ns_O = oxygen.density(lamda);  ns_H = hydrogen.density(lamda); ns_He = helium.density(lamda);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);  wc_He = helium.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda);
        p_mag = sqrt(ppar*ppar+pper*pper);
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix); //<get> efficiency comperable to accessing a member of a struct
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar,pper,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
        vres_f(kz,w_h,alpha,vresz,Eres); //Called only once in first step...
        //Check print parameters
        //std::cout<<"\n"<< R1<<"\n"<< R2<<"\n"<< w1<<"\n"<< w2<<"\n"<< beta<<"\n"<< gama<<"\n"<< mu<<"\n"<< Bywc<<"\n"<< Bzwc<<"\n"<< Bwc<<"\n"<< S<<"\n"<< D<<"\n"<< P<<"\n"<< R<<"\n"<< L<<"\n"<< kappa<<"\n"<< kx<<"\n"<< kz<<"\n"<< w_h<<"\n"<< dwh_ds<<"\n"<< gama;
        
        //RK step-1//#################################################################################################################################################################################################################################################################################################
        slopes(k1, l1, m1, n1, o1, p1, q1, ppar, pper, lamda, eta, alpha, aeq, p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn " << n1<< "\no1 " << o1 << "\np1 " << p1 << "\nq1 " << q1 <<"\n";	
        
        Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o1);
        ns_e = electron.density(lamda+0.5*(Constants::h)*o1); ns_O = oxygen.density(lamda+0.5*(Constants::h)*o1);  ns_H = hydrogen.density(lamda+0.5*(Constants::h)*o1); ns_He = helium.density(lamda+0.5*(Constants::h)*o1);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
        dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o1);
        p_mag = sqrt((ppar+0.5*(Constants::h)*l1)*(ppar+0.5*(Constants::h)*l1)+(pper+0.5*(Constants::h)*m1)*(pper+0.5*(Constants::h)*m1));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar+0.5*(Constants::h)*l1,pper+0.5*(Constants::h)*m1,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
        //RK step-2//#################################################################################################################################################################################################################################################################################################
        slopes(k2, l2, m2, n2, o2, p2, q2, ppar+(0.5*l1*Constants::h), pper+(0.5*m1*Constants::h), lamda+(0.5*o1*Constants::h), eta+(0.5*n1*Constants::h), alpha+(0.5*p1*Constants::h), aeq+(0.5*q1*Constants::h), p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\n";
        

        Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o2);
        ns_e = electron.density(lamda+0.5*(Constants::h)*o2); ns_O = oxygen.density(lamda+0.5*(Constants::h)*o2);  ns_H = hydrogen.density(lamda+0.5*(Constants::h)*o2); ns_He = helium.density(lamda+0.5*(Constants::h)*o2);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
        dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o2);
        p_mag = sqrt((ppar+0.5*(Constants::h)*l2)*(ppar+0.5*(Constants::h)*l2)+(pper+0.5*(Constants::h)*m2)*(pper+0.5*(Constants::h)*m2));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar+0.5*(Constants::h)*l2,pper+0.5*(Constants::h)*m2,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);  
        //RK step-3//#################################################################################################################################################################################################################################################################################################
        slopes(k3, l3, m3, n3, o3, p3, q3, ppar+(0.5*l2*Constants::h), pper+(0.5*m2*Constants::h), lamda+(0.5*o2*Constants::h), eta+(0.5*n2*Constants::h), alpha+(0.5*p2*Constants::h), aeq+(0.5*q2*Constants::h), p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 << "\n";


        Bmag=Bmag_dipole(lamda+(Constants::h)*o3);
        ns_e = electron.density(lamda+(Constants::h)*o3); ns_O = oxygen.density(lamda+(Constants::h)*o3);  ns_H = hydrogen.density(lamda+(Constants::h)*o3); ns_He = helium.density(lamda+(Constants::h)*o3);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda+(Constants::h)*o3);
        p_mag = sqrt((ppar+(Constants::h)*l3)*(ppar+(Constants::h)*l3)+(pper+(Constants::h)*m3)*(pper+(Constants::h)*m3));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar+(Constants::h)*l3,pper+(Constants::h)*m3,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
        //RK step-4//#################################################################################################################################################################################################################################################################################################																								
        slopes(k4, l4, m4, n4, o4, p4, q4, ppar+(l3*Constants::h), pper+(m3*Constants::h), lamda+(o3*Constants::h), eta+(n3*Constants::h), alpha+(p3*Constants::h), aeq+(q3*Constants::h), p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\n";

        
        //Approximate new lamda
        new_lamda = lamda + ((Constants::h)/6)*(o1+2*o2+2*o3+o4);
        
        //Check if NAN to break this particle. Why NAN ?
        if(std::isnan(new_lamda))
        {
            //std::cout<<"\nParticle "<<p<<" breaks";
            //std::cout<<"\n"<< alpha << " " << zeta << " " << ppar<< " " << pper<< " " << eta << " " <<lamda<< " " <<aeq ;
            //std::cout<<"\n" << "ns_He " << ns_He << "\nwc_O " <<wc_O << "\nwc_H " << wc_H << "\nwc_He " << wc_He << "\nwps_e " <<wps_e<< "\nwps_O " <<wps_O << "\nwps_H " << wps_H << "\nwps_He " << wps_He << "\nlamda " <<"\nBwc " << Bwc << "\nEwc "<< Ewc << "\nL " << L << "\nS " <<S<< "\nD " << D << "\nP " << P << "\nR " <<R << "\nmu " << mu << "\nkappa " << kappa<< "\nkx " << kx << "\nkz " <<kz << "\n" << "R1 " << R1 << "\nR2 " << R2 << "\nw1 " << w1 << "\nw2 " << w2 << "\ngama " << gama << "\nbeta " << beta << "Eres " << Eres<< "\nvresz " << vresz << "\nwtau_sq " << wtau_sq << "\nmu_adiabatic" << M_adiabatic;
            break;
        }

        #pragma omp critical //Only one processor can write at a time. There is a chance 2 processors writing in the same spot.
        {                    //This slows down the parallel process, introduces bad scalling 8+ cores. Detecting first and storing in the end demands more memory per process.
            //Check crossing. First estimate new latitude. 
            if( ODPT.crossing(new_lamda*Constants::R2D, lamda*Constants::R2D, Constants::L_shell) )	 
            {										
                //std::cout<<"\nParticle "<< p <<" at: "<<new_lamda*Constants::R2D<< " is about to cross the satellite, at: "<< time << " simulation seconds\n";
                //Store its state(it's before crossing the satellite!).
                ODPT.store( p, lamda, alpha, aeq, time);  			        	
            }
        }

        //Now approximate all values of Runge Kutta's block.
        lamda   =  new_lamda;
        zeta    =  zeta   +  (Constants::h/6)*(k1+2*k2+2*k3+k4);
        ppar    =  ppar   +  (Constants::h/6)*(l1+2*l2+2*l3+l4);
        pper    =  pper   +  (Constants::h/6)*(m1+2*m2+2*m3+m4);
        eta     =  eta    +  (Constants::h/6)*(n1+2*n2+2*n3+n4);
        alpha   =  alpha  +  (Constants::h/6)*(p1+2*p2+2*p3+p4);
        aeq     =  aeq    +  (Constants::h/6)*(q1+2*q2+2*q3+q4);
        //deta_dt =            (Constants::h/6)*(n1+2*n2+2*n3+n4);
        //upar    =  ppar   /  (Constants::m_e*gama);
        //uper    =  pper   /  (Constants::m_e*gama);

        //p_mag = sqrt((ppar*ppar)+(pper*pper));
        //gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //Ekin = ((gama-1)*Constants::m_e*Constants::c*Constants::c)*6.2415e15; 
        //B_lam    =  Bmag_dipole(lamda);    
        //M_adiabatic = (pper*pper)/(2*Constants::m_e*B_lam); 

        //Go to next timestep
        time  = time + Constants::h; 
        
		//To save states:
		//single.save_state(aeq,alpha,lamda,deta_dt,time);

        i++;  
        //std::cout<<"\n\nalpha "<<alpha<<"\nzeta "<< zeta << "\nppar "<< ppar<< "\npper " << pper<< "\neta " << eta << "\nlamda " <<lamda<< "\naeq " <<aeq ;

        //Stop at equator
        //if(eql_dstr[p].lamda.at(i)>0) {	
        //	break;}	
    }

}