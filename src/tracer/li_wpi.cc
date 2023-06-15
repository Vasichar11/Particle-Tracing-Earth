#include "li_wpi.h"


void li_wpi(const int64_t Nsteps_wpi, const int p, std::vector <real> &lat_int, const std::vector <real> &kx_ray, const std::vector <real> &kz_ray, const std::vector <real> &kappa_ray, const std::vector <real> &Bzw, const std::vector <real> &Ezw, const std::vector <real> &Bw_ray, const std::vector <real> &w1, const std::vector <real> &w2, const std::vector <real> &R1, const std::vector <real> &R2, Particles &single, Telescope &ODPT)
{
//---------------------------------------------------- READ RAY HDF5 ----------------------------------------------------//
//Better read out of function and pass vectors as const reference arguments. That way we don't access disc all the time(?)
//---------------------------------------------------- ASSIGN OBJECT VALUES ----------------------------------------------------//
    std::cout.precision(8);                //Output 16 decimal precise
	std::cout<<std::scientific;		        //For e notation representation

    //Assign last particle states from nowpi simulation.
    real lamda      =  single.lamda00;
    real ppar       =  single.ppar00; 
    real pper       =  single.pper00; 
    real alpha      =  single.alpha00; 
    real aeq        =  single.aeq00; 
    real time       =  single.time00;
    real eta        =  single.eta00; 
    real Ekin       =  single.Ekin00; 
    real deta_dt;
    //real zeta       =  single.zeta00; 
    //real upar       =  single.upar00; 
    //real uper       =  single.uper00;
    //std::cout<<"\n\ntime " << time << "\nalpha "<<alpha*Constants::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlamda " <<lamda*Constants::R2D<< "\naeq "<<aeq*Constants::R2D <<"\neta "<<eta*Constants::R2D ;
//------------------------------------------------- LOOP DECLARATIONS -------------------------------------------------//
    int index;                              //To find minimum difference between latitudes
    int i=0;
    real min_lat,max_lat;
    real p_mag,gama,w_h,dwh_ds,kz,Fpar,Fper,Ftheta;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;
    real new_lamda, new_aeq, new_ppar, new_Ekin;
    real min_detadt, mindetadt_time,  max_dEkin, max_dPA, maxdPA_time, maxEkin_time;
    max_dEkin = max_dPA = 0;
    maxEkin_time = maxdPA_time = 0;
    min_detadt = 1; //What would be the minimum detadt value that gives resonance? Only 0?
    mindetadt_time = 0;
//----------------------------------------------------- WPI -----------------------------------------------------------//

    while(i<Nsteps_wpi)          
    {   
        //Take the data from interpolation and return (pulse_dur) of them, moved by "i" each iteration.
        min_lat=*min_element(lat_int.cbegin() + i, lat_int.cbegin() + Constants::puls_dur + i);  //Minimum lat of wave(in pulse duration).
        max_lat=*max_element(lat_int.cbegin() + i, lat_int.cbegin() + Constants::puls_dur + i);  //Max lat of wave.
        //std::cout<<"\n" << lamda*Constants::R2D << " " << min_lat << "" << max_lat << " "<< lat_int.at(i) << " " << lat_int.back(); 
        //std::cout <<"\n\np_mag "<<p_mag<<"\ndwh_ds "<< dwh_ds<<"\nkz "<< kz  <<"\nFpar "<< Fpar <<"\nFper "<< Fper <<"\nFtheta "<< Ftheta <<"\n";  

        f_always(p_mag, gama, w_h, dwh_ds, lamda, ppar, pper); 
        kz = Ftheta = Fpar = Fper = q1 = 0; 
        index = is_in_packet(min_lat, max_lat, lamda, i, lat_int);  //is_in_packet returns "if and where" WPI happens. 
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, q1, kz, pper, ppar, eta, aeq, alpha, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes(k1,  l1,  m1,  n1,  o1,  p1, ppar, pper, alpha, lamda, eta, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz);
        //std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn1 " << n1<< "\no1 " << o1 << "\np1 " << p1 <<"\nq1 " << q1 <<"\n";

        f_always(p_mag, gama, w_h, dwh_ds, lamda+0.5*Constants::h*o1, ppar+0.5*Constants::h*l1, pper+0.5*Constants::h*m1);
        kz = Ftheta = Fpar = Fper = q2 = 0; 
        index = is_in_packet(min_lat, max_lat, lamda+0.5*Constants::h*o1, i, lat_int);
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, q2, kz, pper+0.5*Constants::h*m1, ppar+0.5*Constants::h*l1, eta+0.5*Constants::h*n1, aeq+0.5*Constants::h*q1, alpha+0.5*Constants::h*p1, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes( k2, l2, m2, n2, o2, p2, ppar + 0.5*Constants::h*l1, pper + 0.5*Constants::h*m1, alpha + 0.5*Constants::h*p1, lamda + 0.5*Constants::h*o1, eta + 0.5*Constants::h*n1, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz ); 
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\nq2 "<< q2 <<"\n";

        f_always(p_mag, gama, w_h, dwh_ds, lamda+0.5*Constants::h*o2, ppar+0.5*Constants::h*l2, pper+0.5*Constants::h*m2);
        kz = Ftheta = Fpar = Fper = q3 = 0; 
        index = is_in_packet(min_lat, max_lat, lamda+0.5*Constants::h*o2, i, lat_int);
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, q3, kz, pper+0.5*Constants::h*m2, ppar+0.5*Constants::h*l2, eta+0.5*Constants::h*n2, aeq+0.5*Constants::h*q2, alpha+0.5*Constants::h*p2, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes( k3, l3, m3, n3, o3, p3, ppar + 0.5*Constants::h*l2, pper + 0.5*Constants::h*m2, alpha + 0.5*Constants::h*p2, lamda + 0.5*Constants::h*o2, eta + 0.5*Constants::h*n2, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz );   
        //std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 << "\nq3 "<< q3 <<"\n";

        f_always(p_mag, gama, w_h, dwh_ds, lamda+Constants::h*o3, ppar+Constants::h*l3, pper+Constants::h*m3);
        kz = Ftheta = Fpar = Fper = q4 = 0; 
        index = is_in_packet(min_lat, max_lat, lamda+Constants::h*o3, i, lat_int);
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, q4, kz, pper+Constants::h*m3, ppar+Constants::h*l3, eta+Constants::h*n3, aeq+Constants::h*q3, alpha+Constants::h*p3, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes( k4, l4, m4, n4, o4, p4, ppar + Constants::h*l3, pper + Constants::h*m3, alpha + Constants::h*p3, lamda + Constants::h*o3, eta + Constants::h*n3, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz ); 
        //std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn4 " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\nq4 " << q4 <<"\n";           





        //Runge kutta 4 first estimations:
        //First make the next step increments(lamda,aeq,ppar) that are needed to characterize crossing-validity-trapping. 
        new_lamda = lamda + (Constants::h/6)*(o1+2*o2+2*o3+o4);
        new_aeq   = aeq   + (Constants::h/6)*(q1+2*q2+2*q3+q4);
        new_ppar  = ppar  + (Constants::h/6)*(l1+2*l2+2*l3+l4);
        new_Ekin   = ((gama-1)*Constants::m_e*Constants::c*Constants::c)*6.2415e15;
        //Calculate Energy(???)
        //p_mag  =  sqrt(ppar*ppar + pper*pper);
        //(???)Ekin   = std::sqrt( (p_mag*p_mag*Constants::c*Constants::c) + pow((Constants::c*Constants::c*Constants::m_e),2) ) / 1.602176487E-16;
        


        //Check Validity:
        if(std::isnan(new_lamda*new_aeq*new_ppar))
        {
            single.nan = true;
            single.nan_state(p); //Save the initial state (after noWPI) for the particle and the time that P.A turned out negative.
            std::cout<<"\nParticleV "<<p<<" nan";
            break; 
        }
        //Check Negative P.A:
        if(alpha<0 || aeq<0)
        {
            single.negative = true;
            single.negative_state(p); //Save the initial state (after noWPI) for the particle and the time that P.A turned out negative.
            std::cout<<"\nParticleN "<<p<<" negative p.a";
            break;
        }
        //Check higher than 180 P.A:
        if(alpha>M_PI)
        {
            single.high = true;
            single.high_state(p); //Save the initial state (after noWPI) for the particle and the time that P.A turned out higher than 180.
            std::cout<<"\nParticleH "<<p<<" above 180 p.a";
            break;
        }
        //Check Trapping:
        if( (0<new_aeq && new_aeq<Constants::alpha_lc) || (new_aeq>M_PI-Constants::alpha_lc && new_aeq<M_PI) ) //True if P.A is less than the loss cone angle(for southward particles too).If particle's equator P.A is less than the loss cone angle for this L_shell, then particle is not trapped. hm=100km.
        {
            single.trapped = false;
        }
        //Check Precipitation:
        if(!single.trapped && (ppar*new_ppar<0) ) //Would bounce if ppar is about to change sign.
        {
            single.escaping_state(p, lamda, aeq, alpha, time);
            single.escaped = true;
            std::cout<<"\n\nParticleE "<<p<<" escaped with aeq " <<aeq*Constants::R2D<< " at time " << time ;
            break;
        }

        //Critical Region to push back values in shared memory ODPT object:
        #pragma omp critical //Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot.
        {                    //This slows down the parallel process, introduces bad scalling 8+ cores. Detecting first and storing in the end demands more memory per process.
            //Check Crossing:
            if( ODPT.crossing(new_lamda*Constants::R2D, lamda*Constants::R2D, Constants::L_shell) )	
            {									
                ODPT.store( p, lamda, aeq, alpha, time); //Store its state(it's before crossing the satellite!).		        	
                //std::cout<<"\nParticle "<< p <<" at: "<<new_lamda*Constants::R2D<< " is about to cross the satellite, at: "<< time << " simulation seconds\n";
            }
        }


        //Save state where particle has the maximum(abs) dEkin or dPA
        if(std::abs(new_Ekin - Ekin)>max_dEkin)
        {
            max_dEkin   = new_Ekin - Ekin; 
            maxEkin_time = time;
        }
        if(std::abs(new_aeq - aeq)>max_dPA)
        {
            max_dPA   = new_aeq - aeq*Constants::R2D; 
            maxdPA_time = time;
        }

        //Runge kutta 4 estimations:
        lamda   = new_lamda;
        aeq     = new_aeq;
        ppar    = new_ppar;
        Ekin    = new_Ekin;
        //Rest increments for the next step:
        pper   +=  (Constants::h/6)*(m1+2*m2+2*m3+m4);
        eta    +=  (Constants::h/6)*(n1+2*n2+2*n3+n4);
        deta_dt =  (n1+2*n2+2*n3+n4)/6;
        alpha  +=  (Constants::h/6)*(p1+2*p2+2*p3+p4);
        time  = time + Constants::h; 


        if(std::abs(deta_dt)<min_detadt)
        {
            min_detadt     = std::abs(deta_dt);
            mindetadt_time = time;
        }


        if(lamda*Constants::R2D>=0) //Loop until they reach the equator
        {
            break;
        }

        i++;  
        //std::cout<<"\n\ntime " << time << "\nalpha "<<alpha*Constants::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlamda " <<lamda*Constants::R2D<< "\naeq "<<aeq*Constants::R2D <<"\neta "<<eta*Constants::R2D ;
    
    }
    
    
    
    //Save state when we have deta_dt close to zero, when we have maximum energy diff, and when we have maximum P.A diff(could be the same state?).
    single.save_state( p, max_dEkin, maxEkin_time, max_dPA, maxdPA_time);


}
