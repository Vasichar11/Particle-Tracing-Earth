#include "no_wpi.h"

// Adiabatic motion.
void no_wpi(const int64_t Nsteps_nowpi, int p, Particles &single, Telescope &ODPT)
{
	// std::cout.precision(64);			// Output precision.
	// std::cout<<std::scientific;		// For e notation representation

    // Assign first particle states.
    real latitude = single.latitude0;
    real ppar = single.ppar0;
    real pper = single.pper0;
    real alpha = single.alpha0;
    real aeq = single.aeq0;
    real time = single.time0;
    real Ekin = single.Ekin0;
    // real zeta = single.zeta0;
    // real upar = single.upar0;
    // real uper = single.uper0;
    // std::cout<<"\n\ntime " << time << "\nalpha "<<alpha*Universal::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlatitude " <<latitude*Universal::R2D<< "\naeq "<<aeq*Universal::R2D;

    // Declare function's variables. Once for each particle. When parallel, declare xcore times?
    real new_latitude, new_ppar, new_aeq, new_Ekin;
    real w_h, dwh_ds, Bmag, p_mag, gama;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,o1,o2,o3,o4,p1,p2,p3,p4;
    // Objects for each specie.
    Species electron(Universal::m_e,  Universal::q_e, 1); 
    Species oxygen  (Universal::m_O,  Universal::q_i, 0.006); 
    Species hydrogen(Universal::m_H,  Universal::q_i, 0.94); 
    Species helium  (Universal::m_He, Universal::q_i, 0.054);

    int i=0;
    
    while(i<Nsteps_nowpi) 
    {

        Bmag=Bmag_dipole(latitude);   
        w_h = electron.wc(Bmag); // Cyclotron frequency.
        dwh_ds=dwh_dsf(w_h,latitude);
        p_mag = sqrt(ppar*ppar+pper*pper);
        gama = sqrt((p_mag*p_mag*Universal::c*Universal::c)+(Universal::m_e*Universal::m_e*Universal::c*Universal::c*Universal::c*Universal::c))/(Universal::m_e*Universal::c*Universal::c);
        // RK step-1//#################################################################################################################################################################################################################################################################################################
        slopes(k1, l1, m1, o1, p1, ppar, pper, latitude, w_h, dwh_ds, gama);
        // std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 <<"\no1 " << o1 << "\np1 " << p1 <<"\n";	
        

        Bmag=Bmag_dipole(latitude+0.5*(Simulation::h)*o1);
        w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,latitude+0.5*(Simulation::h)*o1);
        p_mag = sqrt((ppar+0.5*(Simulation::h)*l1)*(ppar+0.5*(Simulation::h)*l1)+(pper+0.5*(Simulation::h)*m1)*(pper+0.5*(Simulation::h)*m1));
        gama = sqrt((p_mag*p_mag*Universal::c*Universal::c)+(Universal::m_e*Universal::m_e*Universal::c*Universal::c*Universal::c*Universal::c))/(Universal::m_e*Universal::c*Universal::c);
        // RK step-2//#################################################################################################################################################################################################################################################################################################
        slopes(k2, l2, m2, o2, p2, ppar+(0.5*l1*Simulation::h), pper+(0.5*m1*Simulation::h), latitude+(0.5*o1*Simulation::h),w_h, dwh_ds, gama);
        // std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 <<"\no2 " << o2 << "\np2 " << p2 << "\n";
        

        Bmag=Bmag_dipole(latitude+0.5*(Simulation::h)*o2);
        w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,latitude+0.5*(Simulation::h)*o2);
        p_mag = sqrt((ppar+0.5*(Simulation::h)*l2)*(ppar+0.5*(Simulation::h)*l2)+(pper+0.5*(Simulation::h)*m2)*(pper+0.5*(Simulation::h)*m2));
        gama = sqrt((p_mag*p_mag*Universal::c*Universal::c)+(Universal::m_e*Universal::m_e*Universal::c*Universal::c*Universal::c*Universal::c))/(Universal::m_e*Universal::c*Universal::c);
        // RK step-3//#################################################################################################################################################################################################################################################################################################
        slopes(k3, l3, m3, o3, p3, ppar+(0.5*l2*Simulation::h), pper+(0.5*m2*Simulation::h), latitude+(0.5*o2*Simulation::h),w_h, dwh_ds, gama);
        // std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3  << "\no3 " << o3 << "\np3 " << p3 << "\n";


        Bmag=Bmag_dipole(latitude+(Simulation::h)*o3);
        w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,latitude+(Simulation::h)*o3);
        p_mag = sqrt((ppar+(Simulation::h)*l3)*(ppar+(Simulation::h)*l3)+(pper+(Simulation::h)*m3)*(pper+(Simulation::h)*m3));
        gama = sqrt((p_mag*p_mag*Universal::c*Universal::c)+(Universal::m_e*Universal::m_e*Universal::c*Universal::c*Universal::c*Universal::c))/(Universal::m_e*Universal::c*Universal::c);
        // RK step-4//#################################################################################################################################################################################################################################################################################################																								
        slopes(k4, l4, m4, o4, p4, ppar+(l3*Simulation::h), pper+(m3*Simulation::h), latitude+(o3*Simulation::h), w_h, dwh_ds, gama);
        // std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 <<"\no4 " << o4 << "\np4 " << p4 << "\n";



        // Runge kutta 4 first estimations:
        // First make the next step increments(latitude,aeq,ppar) that are needed to characterize crossing-validity-trapping. 
        new_latitude = latitude + (Simulation::h/6)*(o1+2*o2+2*o3+o4);
        // !!NOTE!! Adiabatic motion. Aeq stays the same or changes between two values(?) // Is this ok?
        
        int k;
        if(ppar<0) k=1;
        else k=0;
        new_aeq = pow(-1,k) * asin(sin(alpha)*sqrt(Bmag_dipole(0)/Bmag_dipole(new_latitude))) + k*M_PI; 
        //We should have new_aeq   = aeq; and new_Ekin  = Ekin; for adiabatic motion. In practice, aeq changes between two values. One for northward motion and one for southward while on equator. (θ and 180-θ)
        new_Ekin   = ((gama-1)*Universal::m_e*Universal::c*Universal::c)*6.2415e15;
        new_ppar  = ppar  + (Simulation::h/6)*(l1+2*l2+2*l3+l4);


        // Check Validity:
        if(std::isnan(new_latitude*new_aeq*new_ppar))
        {
            single.nan = true;
            single.nan_state(p); //Save the initial state (after noWPI) for the particle and the time that P.A turned out negative.
            std::cout<<"\nParticle(V) "<<p<<" nan";
            break; 
        }
        // Check Negative P.A:
        if(alpha<0 || aeq<0)
        {
            single.negative = true;
            single.negative_state(p); //Save the initial state (after noWPI) for the particle and the time that P.A turned out negative.
            std::cout<<"\nParticle(N) "<<p<<" negative p.a";
            break;
        }
        // Check higher than 180 P.A:
        if(alpha>M_PI)
        {
            single.high = true;
            single.high_state(p); //Save the initial state (after noWPI) for the particle and the time that P.A turned out higher than 180.
            std::cout<<"\nParticle(H) "<<p<<" above 180 p.a";
            break;
        }
        // Check Trapping:
        if( (0<new_aeq && new_aeq<Simulation::alpha_lc) || (new_aeq>M_PI-Simulation::alpha_lc && new_aeq<M_PI) ) //True if P.A is less than the loss cone angle(for southward particles too).If particle's equator P.A is less than the loss cone angle for this L_shell, then particle is not trapped. hm=100km.
        {
            single.trapped = false;
        }
        // Check Precipitation:
        if(!single.trapped && (ppar*new_ppar<0) ) //Would bounce if ppar is about to change sign.
        {
            single.escaping_state(p, latitude, aeq, alpha, time);
            single.escaped = true;
            std::cout<<"\n\nParticle(E) "<<p<<" escaped with aeq " <<aeq*Universal::R2D<< " at time " << time ;
            break;
        }

        // Critical Region to push back values in shared memory ODPT object:

        // Check Crossing:
        if( ODPT.crossing(new_latitude*Universal::R2D, latitude*Universal::R2D, Distribution::L_shell) )	
        {									
            ODPT.store( p, latitude, aeq, alpha, time); // Store its state(it's before crossing the satellite!).		        	
            // std::cout<<"\nParticle "<< p <<" at: "<<new_latitude*Universal::R2D<< " is about to cross the satellite, at: "<< time << " simulation seconds\n";
        }

        
        
        // Runge kutta 4 estimations:
        latitude = new_latitude;
        aeq = new_aeq;
        ppar = new_ppar;
        Ekin = new_Ekin;
        // Rest increments for the next step:
        alpha += (Simulation::h/6)*(p1+2*p2+2*p3+p4);
        pper += (Simulation::h/6)*(m1+2*m2+2*m3+m4);
        
        time = time + Simulation::h; 
        i++;  
       
        // std::cout<<"\n\ntime " << time <<" latitude " << latitude*Universal::R2D<< " Ekin " << Ekin << "\nalpha "<<alpha*Universal::R2D << " aeq0 " << single.aeq0 * Universal::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlatitude " <<latitude*Universal::R2D<< "\naeq "<<aeq*Universal::R2D ;

    }
    

    // Save last state to return values and continue the simulation with wave (if needed). 
    single.latitude_end = latitude;
    single.ppar_end = ppar;
    single.pper_end = pper;
    single.alpha_end = alpha;
    single.aeq_end = aeq;
    single.time_end = time;
    single.Ekin_end = Ekin;
    single.eta_end = Eta_dstr::value; // Particle "gyrophase" initialization 
    // single.zeta_end = zeta;
    // single.upar_end = upar;
    // single.uper_end = uper;
}




