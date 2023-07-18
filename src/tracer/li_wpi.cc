#include "li_wpi.h"


void li_wpi(const int64_t Nsteps_wpi, const int p, std::vector <real> &lat_int, const std::vector <real> &kx_ray, const std::vector <real> &kz_ray, const std::vector <real> &kappa_ray, const std::vector <real> &Bzw, const std::vector <real> &Ezw, const std::vector <real> &Bw_ray, const std::vector <real> &w1, const std::vector <real> &w2, const std::vector <real> &R1, const std::vector <real> &R2, Particles &single, Telescope &ODPT)
{
//---------------------------------------------------- ASSIGN OBJECT VALUES ----------------------------------------------------//
    //std::cout.precision(8); //Output decimal precision
	//std::cout<<std::scientific; //For e notation representation

    //Assign last particle states from nowpi simulation.
    real latitude = single.latitude_end;
    real ppar = single.ppar_end; 
    real pper = single.pper_end; 
    real alpha = single.alpha_end; 
    real aeq = single.aeq_end; 
    real time = single.time_end;
    real eta = single.eta_end; 
    real Ekin = single.Ekin_end; 
    real deta_dt;
    //real zeta = single.zeta_end; 
    //real upar = single.upar_end; 
    //real uper = single.uper_end;
    //std::cout<<"\n\ntime " << time << "\nalpha "<<alpha*Universal::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlatitude " <<latitude*Universal::R2D<< "\naeq "<<aeq*Universal::R2D <<"\neta "<<eta*Universal::R2D ;
//------------------------------------------------- LOOP DECLARATIONS -------------------------------------------------//
    int index; //To find minimum difference between latitudes
    real min_lat , max_lat;
    real p_mag, gama, w_h, dwh_ds, kz, Fpar, Fper, Ftheta;
    real k1, k2, k3, k4; //Runge kutta variables for the slopes
    real l1, l2, l3, l4;
    real m1, m2, m3, m4;
    real n1, n2, n3, n4;
    real o1, o2, o3, o4;
    real p1, p2, p3, p4;
    real q1, q2, q3, q4;
    real new_latitude, new_aeq, new_ppar, new_Ekin;
    real min_detadt, mindetadt_time,  max_dEkin, max_dPA, maxdPA_time, maxEkin_time;
    max_dEkin = max_dPA = 0; // Just a small value
    maxEkin_time = maxdPA_time = mindetadt_time = 0; // Just a small value
    min_detadt = 100; // Just a large value 
//----------------------------------------------------- WPI -----------------------------------------------------------//
    for(int i=0; i<Nsteps_wpi; i++) {
        
        // Take the data from interpolation and return (pulse_dur) of them, moved by "i" each iteration
        min_lat=*min_element(lat_int.cbegin() + i, lat_int.cbegin() + Simulation::puls_dur + i);  //Minimum lat of wave(in pulse duration).
        max_lat=*max_element(lat_int.cbegin() + i, lat_int.cbegin() + Simulation::puls_dur + i);  //Max lat of wave.
        //std::cout<<"\n" << latitude*Universal::R2D << " " << min_lat << "" << max_lat << " "<< lat_int.at(i) << " " << lat_int.back(); 
        //std::cout <<"\n\np_mag "<<p_mag<<"\ndwh_ds "<< dwh_ds<<"\nkz "<< kz  <<"\nFpar "<< Fpar <<"\nFper "<< Fper <<"\nFtheta "<< Ftheta <<"\n";  

        f_always(p_mag, gama, w_h, dwh_ds, latitude, ppar, pper); 
        kz = Ftheta = Fpar = Fper = q1 = 0; 
        index = is_in_packet(min_lat, max_lat, latitude, i, lat_int);  //is_in_packet returns "if and where" WPI happens.
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, pper, ppar, gama, w_h, p_mag, kx_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes(k1, l1, m1, n1, o1, p1, q1, ppar, pper, alpha, latitude, eta, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz, Bw_ray[index], p_mag, aeq, kappa_ray[index]);
        std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn1 " << n1<< "\no1 " << o1 << "\np1 " << p1 <<"\nq1 " << q1 <<"\n";

        f_always(p_mag, gama, w_h, dwh_ds, latitude+0.5*Simulation::h*o1, ppar+0.5*Simulation::h*l1, pper+0.5*Simulation::h*m1);
        kz = Ftheta = Fpar = Fper = q2 = 0; 
        index = is_in_packet(min_lat, max_lat, latitude+0.5*Simulation::h*o1, i, lat_int);
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta,pper+0.5*Simulation::h*m1, ppar+0.5*Simulation::h*l1, gama, w_h, p_mag, kx_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes(k2, l2, m2, n2, o2, p2, q2, ppar + 0.5*Simulation::h*l1, pper + 0.5*Simulation::h*m1, alpha + 0.5*Simulation::h*p1, latitude + 0.5*Simulation::h*o1, eta + 0.5*Simulation::h*n1, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz , Bw_ray[index], p_mag, aeq, kappa_ray[index]); 
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\nq2 "<< q2 <<"\n";

        f_always(p_mag, gama, w_h, dwh_ds, latitude+0.5*Simulation::h*o2, ppar+0.5*Simulation::h*l2, pper+0.5*Simulation::h*m2);
        kz = Ftheta = Fpar = Fper = q3 = 0; 
        index = is_in_packet(min_lat, max_lat, latitude+0.5*Simulation::h*o2, i, lat_int);
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, pper+0.5*Simulation::h*m2, ppar+0.5*Simulation::h*l2, gama, w_h, p_mag, kx_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes(k3, l3, m3, n3, o3, p3, q3, ppar + 0.5*Simulation::h*l2, pper + 0.5*Simulation::h*m2, alpha + 0.5*Simulation::h*p2, latitude + 0.5*Simulation::h*o2, eta + 0.5*Simulation::h*n2, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz , Bw_ray[index], p_mag, aeq, kappa_ray[index]);   
        //std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 << "\nq3 "<< q3 <<"\n";

        f_always(p_mag, gama, w_h, dwh_ds, latitude+Simulation::h*o3, ppar+Simulation::h*l3, pper+Simulation::h*m3);
        kz = Ftheta = Fpar = Fper = q4 = 0; 
        index = is_in_packet(min_lat, max_lat, latitude+Simulation::h*o3, i, lat_int);
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, pper+Simulation::h*m3, ppar+Simulation::h*l3, gama, w_h, p_mag, kx_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes(k4, l4, m4, n4, o4, p4, q4, ppar + Simulation::h*l3, pper + Simulation::h*m3, alpha + Simulation::h*p3, latitude + Simulation::h*o3, eta + Simulation::h*n3, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz , Bw_ray[index], p_mag, aeq, kappa_ray[index]); 
        //std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn4 " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\nq4 " << q4 <<"\n";           


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
        if (std::isnan(new_latitude) || std::isnan(new_aeq) || std::isnan(new_ppar)) {
            single.nan = true;
            single.nan_state(p); // Save id of particle
            std::cout << "\nParticleV " << p << " develops NaN state";
            std::cout << "\nq1+2*q2+2*q3+q4 " << q1<<" "<<q2<<" "<< q3<< " " <<q4;
            break; 
        }
        // Check Negative P.A:
        if(alpha<0 || aeq<0) {
            single.negative = true;
            single.negative_state(p); // Save id of particle
            std::cout<<"\nParticleN "<<p<<" negative P.A";
            break;
        }
        // Check higher than 180 P.A:
        if(alpha>M_PI) {
            single.high = true;
            single.high_state(p); // Save id of particle
            std::cout<<"\nParticleH "<<p<<" above 180 P.A";
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
            std::cout<<"\n\nParticleE "<<p<<" escaped with aeq " <<aeq*Universal::R2D<< " at time " << time ;
            break;
        }
        // Critical Region to push back values in shared memory ODPT object
        // Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot
        // This can slow down the parallel process, and will introduce bad scalling 8+ cores, which is sufficient for this simulation
        // Alternative: Detecting first and storing in the end demands more memory per process
        #pragma omp critical
        {                    
            //Check Crossing:
            if( ODPT.crossing(new_latitude*Universal::R2D, latitude*Universal::R2D, Distribution::L_shell) ) {									
                ODPT.store( p, latitude, aeq, alpha, time); //Store its state (just before crossing the satellite)	        	
                //std::cout<<"\nParticle "<< p <<" at: "<<new_latitude*Universal::R2D<< " is about to cross the satellite, at: "<< time << " simulation seconds\n";
            }
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
        deta_dt = (n1+2*n2+2*n3+n4)/6;
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

        std::cout<<"\n\ntime " << time << "\nalpha "<<alpha*Universal::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlatitude " <<latitude*Universal::R2D<< "\naeq "<<aeq*Universal::R2D <<"\neta "<<eta*Universal::R2D ;
        std::cout<<"\naeq"<<aeq;
    
    }
    
    
    
    // Save state when we have deta_dt close to zero, when we have maximum energy diff, and when we have maximum P.A diff.
    // single.save_state( p, max_dEkin, maxEkin_time, max_dPA, maxdPA_time, min_detadt);


}

