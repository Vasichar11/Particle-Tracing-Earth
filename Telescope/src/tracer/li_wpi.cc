#include "headers/li_wpi.h"


void li_wpi(real p, Particles &single, Telescope &ODPT)
{
//---------------------------------------------------- READ RAY HDF5 ----------------------------------------------------//
    //read_vector() is a function to read HDF5 dataset as vectors. 
    static std::vector <real> lat_int       =   read_vector("lat_int",       "h5files/interpolated_ray.h5");
    static std::vector <real> kx_ray        =   read_vector("kx_ray",        "h5files/interpolated_ray.h5");    
    static std::vector <real> kz_ray        =   read_vector("kz_ray",        "h5files/interpolated_ray.h5");   
    static std::vector <real> kappa_ray     =   read_vector("kappa_ray",     "h5files/interpolated_ray.h5");       
    static std::vector <real> Bzw           =   read_vector("Bzw",           "h5files/interpolated_ray.h5");
    static std::vector <real> Ezw           =   read_vector("Ezw",           "h5files/interpolated_ray.h5");
    static std::vector <real> Bw_ray        =   read_vector("Bw_ray",        "h5files/interpolated_ray.h5");    
    static std::vector <real> w1            =   read_vector("w1",            "h5files/interpolated_ray.h5");
    static std::vector <real> w2            =   read_vector("w2",            "h5files/interpolated_ray.h5");
    static std::vector <real> R1            =   read_vector("R1",            "h5files/interpolated_ray.h5");
    static std::vector <real> R2            =   read_vector("R2",            "h5files/interpolated_ray.h5");
//---------------------------------------------------- ASSIGN OBJECT VALUES ----------------------------------------------------//

    


    //Assign last particle states from nowpi simulation.
    real lamda    =  single.lamda_end;
    real ppar     =  single.ppar_end; 
    real pper     =  single.pper_end; 
    real alpha    =  single.alpha_end; 
    real aeq      =  single.aeq_end; 
    real time     =  single.time_end;
    real eta      =  single.eta_end; 
    //real zeta     =  single.zeta_end; 
    //real upar     =  single.upar_end; 
    //real uper     =  single.uper_end;
//------------------------------------------------- LOOP DECLARATIONS -------------------------------------------------//
    int index;                              //To find minimum difference between latitudes
    int i=0;

    real min_lat,max_lat;
    real p_mag,gama,w_h,dwh_ds,kz,Fpar,Fper,Ftheta;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;
    real new_lamda, new_aeq, new_ppar;
    bool trapped = true;                       //Particles trapped in Earth's magnetic field.

    //std::cout.precision(8);                //Output 16 decimal precise
	//std::cout<<std::scientific;		        //For e notation representation
//----------------------------------------------------- WPI -----------------------------------------------------------//

    while(i<Constants::Nsteps_wpi)          
    {   
        //Take the data from interpolation and return (pulse_dur) of them, moved by "i" each iteration.
        min_lat=*min_element(lat_int.cbegin() + i, lat_int.cbegin() + Constants::puls_dur + i);  //Minimum lat of wave(in pulse duration).
        max_lat=*max_element(lat_int.cbegin() + i, lat_int.cbegin() + Constants::puls_dur + i);  //Max lat of wave.
        //std::cout<<"\n" << lamda*Constants::R2D << " " << min_lat << "" << max_lat << " "<< lat_int.at(i) << " " << lat_int.back(); 
        
        f_always(p_mag, gama, w_h, dwh_ds, lamda, ppar, pper); 
        kz = Ftheta = Fpar = Fper = q1 = 0; 
        index = is_in_packet(min_lat, max_lat, lamda, i, lat_int);  //is_in_packet returns "if and where" WPI happens. 
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, q1, kz, pper, ppar, eta, aeq, alpha, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes(k1,  l1,  m1,  n1,  o1,  p1, ppar, pper, alpha, lamda, eta, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz);
        //std::cout <<"\np_mag "<<p_mag<<"\ndwh_ds "<< dwh_ds<<"\nkz "<< kz  <<"\nFpar "<< Fpar <<"\nFper "<< Fper <<"\nFtheta "<< Ftheta <<"\n";  
        //std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn1 " << n1<< "\no1 " << o1 << "\np1 " << p1 <<"\nq1 " << q1 <<"\n";

        f_always(p_mag, gama, w_h, dwh_ds, lamda+0.5*Constants::h*o1, ppar+0.5*Constants::h*l1, pper+0.5*Constants::h*m1);
        kz = Ftheta = Fpar = Fper = q2 = 0; 
        index = is_in_packet(min_lat, max_lat, lamda+0.5*Constants::h*o1, i, lat_int);
        if(index>=0) { //Do only if there's WPI
            f_packet(Fpar, Fper, Ftheta, q2, kz, pper+0.5*Constants::h*m1, ppar+0.5*Constants::h*l1, eta+0.5*Constants::h*n1, aeq+0.5*Constants::h*q1, alpha+0.5*Constants::h*p1, gama, w_h, p_mag, kx_ray[index], kz_ray[index], kappa_ray[index], Bw_ray[index], Bzw[index], Ezw[index], w1[index], w2[index], R1[index], R2[index]);
        }
        slopes( k2, l2, m2, n2, o2, p2, ppar + 0.5*Constants::h*l1, pper + 0.5*Constants::h*m1, alpha + 0.5*Constants::h*p1, lamda + 0.5*Constants::h*o1, eta + 0.5*Constants::h*n1, Fpar, Fper, Ftheta, gama, w_h, dwh_ds, kz ); 
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\nq2 "<< q2 <<"\n";;
        
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


        //Check Validity:
        new_lamda = lamda + (Constants::h/6)*(o1+2*o2+2*o3+o4); //Approximate new lamda first
        if(std::isnan(new_lamda)) { std::cout<<"\nParticle "<<p<<" breaks"; break; }
        if(alpha<0 || aeq<0)      { std::cout<<"\nParticle "<<p<<" negative p.a"; break; }
        
        //Check Crossing:
        #pragma omp critical //Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot.
        {                    //This slows down the parallel process, introduces bad scalling 8+ cores. Detecting first and storing in the end demands more memory per process.
            if( ODPT.crossing(new_lamda*Constants::R2D, lamda*Constants::R2D, Constants::L_shell) )	 
            {	//Check crossing.								
                //std::cout<<"\nParticle "<< p <<" at: "<<new_lamda*Constants::R2D<< " is about to cross the satellite, at: "<< time << " simulation seconds\n";
                ODPT.store( p, lamda, alpha, aeq, time); //Store its state(it's before crossing the satellite!).		        	
            }
        }
        //Check Trapping:
        new_aeq = aeq + (Constants::h/6)*(q1+2*q2+2*q3+q4);
        if( (0<new_aeq && new_aeq<Constants::alpha_lc) || (new_aeq>M_PI-Constants::alpha_lc && new_aeq<M_PI) ) //True if P.A is less than the loss cone angle(for southward particles too).
        {                                                //If particle's equator P.A is less than the loss cone angle for this L_shell, then particle is not trapped. hm=100km.
            trapped = false;
        }
        //Check Precipitation:
        new_ppar = ppar + (Constants::h/6)*(l1+2*l2+2*l3+l4);
        if(!trapped && (ppar*new_ppar<0) ) //Would bounce if ppar is about to change sign.
        {   
            //To save states of precipitating particles:
            #pragma omp critical //Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot.
            {   
                single.save_state(p, lamda, alpha, aeq, time);
                //std::cout<<"\n\nParticle "<<p<<" escaped with ppar "<<ppar<< " new_ppar would be "<<new_ppar<<" pper " << pper<< " eta " << eta << " lamda " <<lamda*Constants::R2D<< " alpha "<< alpha*Constants::R2D << " aeq " <<aeq*Constants::R2D<< " at time " << time ;
                single.escaped = true;
            }
            break;
        }
        //Next step:
        new_values_RK4(lamda, ppar, pper, eta, alpha, aeq, l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4, o1, o2, o3, o4, p1, p2, p3, p4, q1, q2, q3, q4);
        time  = time + Constants::h; 
        i++;  

        //std::cout<<"\n\nppar "<< ppar<< "\npper " << pper<< "\neta " << eta << "\nlamda " <<lamda*Constants::R2D<< "\nalpha "<< alpha*Constants::R2D << "\naeq " <<aeq*Constants::R2D ;
    }


}

