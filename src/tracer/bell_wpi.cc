#include "bell_wpi.h"

// For wave-particle interaction using Bell equations
// For now, the wave is considered to be everywhere
void bell_wpi(const int64_t Nsteps_wpi, int p, Particles &single, Telescope &ODPT)
{
	//std::cout.precision(64); // Output 16 decimal precise
	//std::cout<<std::scientific;// For e notation representation
    
    real latitude_end = single.latitude_end;
    real ppar_end = single.ppar_end; 
    real pper_end = single.pper_end; 
    real alpha_end = single.alpha_end; 
    real aeq_end = single.aeq_end; 
    real eta_end = single.eta_end; 
    real time_end = single.time_end;
    real Ekin_end = single.Ekin_end; 
    //real zeta_end = single.zeta_end; 
    //real upar_end = single.upar_end; 
    //real uper_end = single.uper_end;
    
    real latitude = latitude_end;
    real ppar = ppar_end; 
    real pper = pper_end; 
    real alpha = alpha_end; 
    real aeq = aeq_end; 
    real eta = eta_end; 
    real time = time_end;
    real Ekin = Ekin_end;
    //real zeta = zeta_end; 
    //real upar = upar_end; 
    //real uper = uper_end;

    //std::cout<<"\n\nParticle "<<p<<" at alpha "<<alpha << "\nppar "<< ppar<< "\npper " << pper<< "\neta " << eta << "\nlatitude " <<latitude<< "\naeq " <<aeq ;

    // Declare function's variables. Once for each particle
    real new_latitude = latitude;
    real new_ppar, new_aeq, new_Ekin;
    real ns_e,w_h, wps_e, ns_O, wc_O, wps_O ,ns_H, wc_H, wps_H, ns_He, wc_He, wps_He;
    real Bmag;
    real k1, k2, k3, k4;
    real l1, l2, l3, l4;
    real m1, m2, m3, m4;
    real n1, n2, n3, n4;
    real o1, o2, o3, o4;
    real p1, p2, p3, p4;
    real q1, q2, q3, q4;
    real l1_old, l2_old, l3_old, l4_old;
    real m1_old, m2_old, m3_old, m4_old;
    real n1_old, n2_old, n3_old, n4_old;
    real o1_old, o2_old, o3_old, o4_old;
    real p1_old, p2_old, p3_old, p4_old;
    real q1_old, q2_old, q3_old, q4_old;
    l1_old=l2_old=l3_old=l4_old=m1_old=m2_old=m3_old=m4_old=n1_old=n2_old=n3_old=n4_old=o1_old=o2_old=o3_old=o4_old=p1_old=p2_old=p3_old=p4_old=q1_old=q2_old=q3_old=q4_old=0;
    real gama, w1, w2, R1, R2, beta, wtau_sq;
    real S, D, P, R, L, mu, dwh_ds, kappa, kx, kz;
    // real vres, Eres;
    real Bxwc, Bzwc, Bywc, Exwc, Eywc, Ezwc, Bwc;
    real p_mag;
    real min_detadt, mindetadt_time,  max_dEkin, max_dPA, maxdPA_time, maxEkin_time;
    max_dEkin = max_dPA = 0; // Just a small value
    maxEkin_time = maxdPA_time = mindetadt_time = 0; // Just a small value
    min_detadt = 100; // Just a large value 

    // Tuples for WPI
    std::tuple<real, real, real, real, real> stix;
    std::tuple<real, real, real, real> disp;

    // Objects for each specie
    Species electron(Universal::m_e,  Universal::q_e, 1); 
    Species oxygen  (Universal::m_O,  Universal::q_i, 0.006); 
    Species hydrogen(Universal::m_H,  Universal::q_i, 0.94); 
    Species helium  (Universal::m_He, Universal::q_i, 0.054);
    


    for(int i=0; i<Nsteps_wpi; i++) {
        
        Bmag=Bmag_dipole(latitude);
        ns_e = electron.density(latitude); ns_O = oxygen.density(latitude);  ns_H = hydrogen.density(latitude); ns_He = helium.density(latitude);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);  wc_He = helium.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,latitude);
        p_mag = sqrt(ppar*ppar+pper*pper);
        gama = sqrt((p_mag*p_mag*Universal::c*Universal::c)+(Universal::m_e*Universal::m_e*Universal::c*Universal::c*Universal::c*Universal::c))/(Universal::m_e*Universal::c*Universal::c);
        // For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix); //<get> efficiency comperable to accessing a member of a struct
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar,pper,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
        //vres_f(kz,w_h,alpha,vresz,Eres); //Called only once in first step...
        //std::cout<<"\nR1 "<< R1<<"\nR2 "<< R2<<"\nw1 "<< w1<<"\nw2 "<< w2<<"\nbeta "<< beta<<"\ngama "<< gama<<"\nmu "<< mu<<"\nBywc "<< Bywc<<"\nBzwc "<< Bzwc<<"\nBwc "<< Bwc<<"\nS "<< S<<"\nD "<< D<<"\nP "<< P<<"\nR "<< R<<"\nL "<< L<<"\nkappa "<< kappa<<"\nkx "<< kx<<"\nkz "<< kz<<"\nw_h "<< w_h<<"\ndwh_ds "<< dwh_ds<<"\ngama "<< gama;
        slopes(k1, l1, m1, n1, o1, p1, q1, ppar, pper, latitude, eta, alpha, aeq, p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn " << n1<< "\no1 " << o1 << "\np1 " << p1 << "\nq1 " << q1 <<"\n";	


        Bmag=Bmag_dipole(latitude+0.5*(Simulation::h)*o1);
        ns_e = electron.density(latitude+0.5*(Simulation::h)*o1); ns_O = oxygen.density(latitude+0.5*(Simulation::h)*o1);  ns_H = hydrogen.density(latitude+0.5*(Simulation::h)*o1); ns_He = helium.density(latitude+0.5*(Simulation::h)*o1);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
        dwh_ds=dwh_dsf(w_h,latitude+0.5*(Simulation::h)*o1);
        p_mag = sqrt((ppar+0.5*(Simulation::h)*l1)*(ppar+0.5*(Simulation::h)*l1)+(pper+0.5*(Simulation::h)*m1)*(pper+0.5*(Simulation::h)*m1));
        gama = sqrt((p_mag*p_mag*Universal::c*Universal::c)+(Universal::m_e*Universal::m_e*Universal::c*Universal::c*Universal::c*Universal::c))/(Universal::m_e*Universal::c*Universal::c);
        // For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar+0.5*(Simulation::h)*l1,pper+0.5*(Simulation::h)*m1,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);
        //std::cout<<"\nR1 "<< R1<<"\nR2 "<< R2<<"\nw1 "<< w1<<"\nw2 "<< w2<<"\nbeta "<< beta<<"\ngama "<< gama<<"\nmu "<< mu<<"\nBywc "<< Bywc<<"\nBzwc "<< Bzwc<<"\nBwc "<< Bwc<<"\nS "<< S<<"\nD "<< D<<"\nP "<< P<<"\nR "<< R<<"\nL "<< L<<"\nkappa "<< kappa<<"\nkx "<< kx<<"\nkz "<< kz<<"\nw_h "<< w_h<<"\ndwh_ds "<< dwh_ds<<"\ngama "<< gama;
        slopes(k2, l2, m2, n2, o2, p2, q2, ppar+(0.5*l1*Simulation::h), pper+(0.5*m1*Simulation::h), latitude+(0.5*o1*Simulation::h), eta+(0.5*n1*Simulation::h), alpha+(0.5*p1*Simulation::h), aeq+(0.5*q1*Simulation::h), p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 <<"\nq2 "<< q2 <<"\n";
        

        Bmag=Bmag_dipole(latitude+0.5*(Simulation::h)*o2);
        ns_e = electron.density(latitude+0.5*(Simulation::h)*o2); ns_O = oxygen.density(latitude+0.5*(Simulation::h)*o2);  ns_H = hydrogen.density(latitude+0.5*(Simulation::h)*o2); ns_He = helium.density(latitude+0.5*(Simulation::h)*o2);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
        dwh_ds=dwh_dsf(w_h,latitude+0.5*(Simulation::h)*o2);
        p_mag = sqrt((ppar+0.5*(Simulation::h)*l2)*(ppar+0.5*(Simulation::h)*l2)+(pper+0.5*(Simulation::h)*m2)*(pper+0.5*(Simulation::h)*m2));
        gama = sqrt((p_mag*p_mag*Universal::c*Universal::c)+(Universal::m_e*Universal::m_e*Universal::c*Universal::c*Universal::c*Universal::c))/(Universal::m_e*Universal::c*Universal::c);
        // For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar+0.5*(Simulation::h)*l2,pper+0.5*(Simulation::h)*m2,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);  
        //std::cout<<"\nR1 "<< R1<<"\nR2 "<< R2<<"\nw1 "<< w1<<"\nw2 "<< w2<<"\nbeta "<< beta<<"\ngama "<< gama<<"\nmu "<< mu<<"\nBywc "<< Bywc<<"\nBzwc "<< Bzwc<<"\nBwc "<< Bwc<<"\nS "<< S<<"\nD "<< D<<"\nP "<< P<<"\nR "<< R<<"\nL "<< L<<"\nkappa "<< kappa<<"\nkx "<< kx<<"\nkz "<< kz<<"\nw_h "<< w_h<<"\ndwh_ds "<< dwh_ds<<"\ngama "<< gama;
        slopes(k3, l3, m3, n3, o3, p3, q3, ppar+(0.5*l2*Simulation::h), pper+(0.5*m2*Simulation::h), latitude+(0.5*o2*Simulation::h), eta+(0.5*n2*Simulation::h), alpha+(0.5*p2*Simulation::h), aeq+(0.5*q2*Simulation::h), p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 <<"\nq3 "<< q3 <<"\n";
        

        Bmag=Bmag_dipole(latitude+Simulation::h*o3);
        ns_e = electron.density(latitude+Simulation::h*o3); ns_O = oxygen.density(latitude+Simulation::h*o3);  ns_H = hydrogen.density(latitude+Simulation::h*o3); ns_He = helium.density(latitude+Simulation::h*o3);
        w_h = electron.wc(Bmag); wc_O = oxygen.wc(Bmag); wc_H = hydrogen.wc(Bmag);   wc_He = helium.wc(Bmag); 
        dwh_ds=dwh_dsf(w_h,latitude+Simulation::h*o3);
        p_mag = sqrt((ppar+Simulation::h*l3)*(ppar+Simulation::h*l3)+(pper+Simulation::h*m3)*(pper+Simulation::h*m3));
        gama = sqrt((p_mag*p_mag*Universal::c*Universal::c)+(Universal::m_e*Universal::m_e*Universal::c*Universal::c*Universal::c*Universal::c))/(Universal::m_e*Universal::c*Universal::c);
        // For interaction
        wps_e = electron.wps(ns_e); wps_O = oxygen.wps(ns_O); wps_H = hydrogen.wps(ns_H); wps_He = helium.wps(ns_He); 
        stix = stix_parameters(w_h, wc_O, wc_H, wc_He, wps_e, wps_O, wps_H, wps_He);
        S=std::get<0>(stix); D=std::get<1>(stix); P=std::get<2>(stix); R=std::get<3>(stix); L=std::get<4>(stix);
        disp = dispersion(S,P,R,L,D);
        mu=std::get<0>(disp); kappa=std::get<1>(disp); kx=std::get<2>(disp); kz=std::get<3>(disp);
        whistlers(p,i,mu,P,D,S,kz, Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc);
        //Ewc = sqrt(Exwc*Exwc + Eywc*Eywc + Ezwc*Ezwc);
        Bwc = sqrt(Bxwc*Bxwc + Bywc*Bywc + Bzwc*Bzwc);
        Bell_params(ppar+Simulation::h*l3,pper+Simulation::h*m3,Bxwc,Bywc,Exwc,Eywc,Ezwc,kz,kx,w_h,gama,w1,w2,wtau_sq,R1,R2,beta);  
        //std::cout<<"\nR1 "<< R1<<"\nR2 "<< R2<<"\nw1 "<< w1<<"\nw2 "<< w2<<"\nbeta "<< beta<<"\ngama "<< gama<<"\nmu "<< mu<<"\nBywc "<< Bywc<<"\nBzwc "<< Bzwc<<"\nBwc "<< Bwc<<"\nS "<< S<<"\nD "<< D<<"\nP "<< P<<"\nR "<< R<<"\nL "<< L<<"\nkappa "<< kappa<<"\nkx "<< kx<<"\nkz "<< kz<<"\nw_h "<< w_h<<"\ndwh_ds "<< dwh_ds<<"\ngama "<< gama;
        slopes(k4, l4, m4, n4, o4, p4, q4, ppar+l3*Simulation::h, pper+m3*Simulation::h, latitude+o3*Simulation::h, eta+n3*Simulation::h, alpha+p3*Simulation::h, aeq+q3*Simulation::h, p_mag, w_h, dwh_ds, gama, kz, kappa, wtau_sq, w1, w2, R1, R2, beta, Bwc);
        //std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\nq4 "<< q4 <<"\n";
       
        
        // Old slope values kept in memory to encounter the NAN case
        // NAN case is probably caused because of precision issues
        // isnan(mu) --> step once more using old slopes to reach to a valid state
        // i.e when NAN, step of h becomes 2*h, 3*h ... until valid
        l1_old = l1 ; l2_old = l2 ; l3_old = l3 ; l4_old = l4 ;  
        m1_old = m1 ; m2_old = m2 ; m3_old = m3 ; m4_old = m4 ;  
        n1_old = n1 ; n2_old = n2 ; n3_old = n3 ; n4_old = n4 ;  
        o1_old = o1 ; o2_old = o2 ; o3_old = o3 ; o4_old = o4 ;  
        p1_old = p1 ; p2_old = p2 ; p3_old = p3 ; p4_old = p4 ;  
        q1_old = q1 ; q2_old = q2 ; q3_old = q3 ; q4_old = q4 ; 

        // Runge kutta 4 first estimations:
        // First make the next step increments(latitude,aeq,ppar) that are needed to characterize crossing-validity-trapping 
        new_latitude = latitude + (Simulation::h/6)*(o1+2*o2+2*o3+o4);
        new_aeq = aeq + (Simulation::h/6)*(q1+2*q2+2*q3+q4);
        new_ppar = ppar + (Simulation::h/6)*(l1+2*l2+2*l3+l4);
        new_Ekin = ((gama-1)*Universal::m_e*Universal::c*Universal::c)*6.2415e15;
        // Calculate Energy
        // p_mag  =  sqrt(ppar*ppar + pper*pper);
        // Ekin   = std::sqrt( (p_mag*p_mag*Universal::c*Universal::c) + pow((Universal::c*Universal::c*Universal::m_e),2) ) / 1.602176487E-16;

        
        // Check Validity:
        if(std::isnan(new_latitude*new_aeq*new_ppar)) {
            single.nan = true;
            single.nan_state(p); // Save id of particle
            std::cout<<"\nParticle(V) "<<p<<" nan";
            break; 
        }
        // Check Negative P.A:
        if(alpha<0 || aeq<0) {
            single.negative = true;
            single.negative_state(p); // Save id of particle
            std::cout<<"\nParticle(N) "<<p<<" negative p.a";
            break;
        }
        // Check higher than 180 P.A:
        if(alpha>M_PI) {
            single.high = true;
            single.high_state(p); // Save id of particle
            std::cout<<"\nParticle(H) "<<p<<" above 180 p.a";
            break;
        }
        // Check Trapping:
        // True if P.A is less than the loss cone angle(for southward particles too)
        // While aeq<loss cone P.A for this L_shell the particle continues oscillating
        if( (0<new_aeq && new_aeq<Simulation::alpha_lc) || (new_aeq>M_PI-Simulation::alpha_lc && new_aeq<M_PI) ) {
            single.trapped = false;
            // If no longer trapped, don't break immediately
            // The particle will follow its trajectory
            // until reaching the moment just before changing its direction of motion, at which point it will precipitate.
        }
        // Check Precipitation:
        // If no longer trapped, and is about to change direction
        if(!single.trapped && (ppar*new_ppar<0) ) {
            single.escaping_state(p, latitude, aeq, alpha, time);
            single.escaped = true;
            std::cout<<"\n\nParticle(E) "<<p<<" escaped with aeq " <<aeq*Universal::R2D<< " at time " << time ;
            break;
        }
        // Critical Region to push back values in shared memory ODPT object
        // Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot
        // This can slow down the parallel process, and will introduce bad scalling 8+ cores, which is sufficient for this simulation
        // Alternative: Detecting first and storing in the end demands more memory per process
        //Check Crossing:
        if( ODPT.crossing(new_latitude*Universal::R2D, latitude*Universal::R2D, Distribution::L_shell) ) {									
            ODPT.store( p, latitude, aeq, alpha, time); //Store its state (just before crossing the satellite)	        	
            //std::cout<<"\nParticle "<< p <<" at: "<<new_latitude*Universal::R2D<< " is about to cross the satellite, at: "<< time << " simulation seconds\n";
        }
        // Save state where particle has the maximum(abs) dEkin or dPA
        // if(std::abs(new_Ekin - Ekin)>max_dEkin) {
        //     max_dEkin = new_Ekin - Ekin; 
        //     maxEkin_time = time;
        // }
        // if(std::abs(new_aeq - aeq)>max_dPA) {
        //     max_dPA = new_aeq - aeq*Universal::R2D; 
        //     maxdPA_time = time;
        // }

        // Runge kutta 4 estimations:
        latitude = new_latitude;
        aeq = new_aeq;
        ppar = new_ppar;
        Ekin = new_Ekin;
        // Rest increments for the next step:
        pper += (Simulation::h/6)*(m1+2*m2+2*m3+m4);
        eta += (Simulation::h/6)*(n1+2*n2+2*n3+n4);
        // deta_dt = (n1+2*n2+2*n3+n4)/6;
        alpha += (Simulation::h/6)*(p1+2*p2+2*p3+p4);
        time = time + Simulation::h; 


        // if(std::abs(deta_dt)<min_detadt) {
        //     min_detadt = std::abs(deta_dt);
        //     mindetadt_time = time;
        // }

        // Loop until they reach the equator
        // if(latitude*Universal::R2D>=0) {
        //     break;
        // }

        //std::cout<<"\n\nalpha "<<alpha*Universal::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlatitude " <<latitude*Universal::R2D<< "\naeq "<<aeq*Universal::R2D;

    }

    // Save state when we have deta_dt close to zero, when we have maximum energy diff, and when we have maximum P.A diff.
    //single.save_state( p, max_dEkin, maxEkin_time, max_dPA, maxdPA_time, min_detadt);

}
