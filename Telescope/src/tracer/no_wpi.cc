#include "headers/no_wpi.h"

//Adiabatic motion.
void no_wpi(int p, Particles &single, Telescope &ODPT)
{
        
	std::cout.precision(64);			//Output 16 decimal precise
	std::cout<<std::scientific;		//For e notation representation

    real lamda    =  single.lamda_init;
    real ppar     =  single.ppar_init; 
    real pper     =  single.pper_init; 
    real alpha    =  single.alpha_init; 
    real aeq      =  single.aeq_init; 
    real time     =  single.time_init;
    //real zeta     =  single.zeta_init; 
    //real upar     =  single.upar_init; 
    //real uper     =  single.uper_init;
    std::cout<<"\n\nalpha "<<alpha*Constants::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlamda " <<lamda*Constants::R2D<< "\naeq "<<aeq*Constants::R2D;

    //Declare function's variables. Once for each particle. When parallel, declare xcore times?
    real new_lamda, new_ppar;
    real w_h, dwh_ds, Bmag, p_mag, gama;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,o1,o2,o3,o4,p1,p2,p3,p4;

    //Objects for each specie.
    Species electron(Constants::m_e,  Constants::q_e, 1); 
    Species oxygen  (Constants::m_O,  Constants::q_i, 0.006); 
    Species hydrogen(Constants::m_H,  Constants::q_i, 0.94); 
    Species helium  (Constants::m_He, Constants::q_i, 0.054);

    int i=0;
    
    while(i<Constants::Nsteps_nowpi) 
    {

        Bmag=Bmag_dipole(lamda);   
        w_h = electron.wc(Bmag); //Cyclotron frequency.
        dwh_ds=dwh_dsf(w_h,lamda);
        p_mag = sqrt(ppar*ppar+pper*pper);
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-1//#################################################################################################################################################################################################################################################################################################
        slopes(k1, l1, m1, o1, p1, ppar, pper, lamda, w_h, dwh_ds, gama);
        //std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 <<"\no1 " << o1 << "\np1 " << p1 <<"\n";	
        

        Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o1);
        w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o1);
        p_mag = sqrt((ppar+0.5*(Constants::h)*l1)*(ppar+0.5*(Constants::h)*l1)+(pper+0.5*(Constants::h)*m1)*(pper+0.5*(Constants::h)*m1));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-2//#################################################################################################################################################################################################################################################################################################
        slopes(k2, l2, m2, o2, p2, ppar+(0.5*l1*Constants::h), pper+(0.5*m1*Constants::h), lamda+(0.5*o1*Constants::h),w_h, dwh_ds, gama);
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 <<"\no2 " << o2 << "\np2 " << p2 << "\n";
        

        Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o2);
        w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o2);
        p_mag = sqrt((ppar+0.5*(Constants::h)*l2)*(ppar+0.5*(Constants::h)*l2)+(pper+0.5*(Constants::h)*m2)*(pper+0.5*(Constants::h)*m2));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-3//#################################################################################################################################################################################################################################################################################################
        slopes(k3, l3, m3, o3, p3, ppar+(0.5*l2*Constants::h), pper+(0.5*m2*Constants::h), lamda+(0.5*o2*Constants::h),w_h, dwh_ds, gama);
        //std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3  << "\no3 " << o3 << "\np3 " << p3 << "\n";


        Bmag=Bmag_dipole(lamda+(Constants::h)*o3);
        w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda+(Constants::h)*o3);
        p_mag = sqrt((ppar+(Constants::h)*l3)*(ppar+(Constants::h)*l3)+(pper+(Constants::h)*m3)*(pper+(Constants::h)*m3));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-4//#################################################################################################################################################################################################################################################################################################																								
        slopes(k4, l4, m4, o4, p4, ppar+(l3*Constants::h), pper+(m3*Constants::h), lamda+(o3*Constants::h), w_h, dwh_ds, gama);
        //std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 <<"\no4 " << o4 << "\np4 " << p4 << "\n";

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
        if( (0<aeq && aeq<Constants::alpha_lc) || (aeq>M_PI-Constants::alpha_lc && aeq<M_PI) ) //True if P.A is less than the loss cone angle(for southward particles too).
        {                                                //If particle's equator P.A is less than the loss cone angle for this L_shell, then particle is not trapped. hm=100km.
            single.trapped = false;
        }

        //Check Precipitation:
        new_ppar = ppar + (Constants::h/6)*(l1+2*l2+2*l3+l4);
        if(!single.trapped && (ppar*new_ppar<0) ) //Would bounce if ppar is about to change sign.
        {   
            //To save states of precipitating particles:
            #pragma omp critical //Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot.
            {   
                single.escaping_state(p, lamda, alpha, aeq, time);
                //std::cout<<"\n\nParticle "<<p<<" escaped with ppar "<<ppar<< " new_ppar would be "<<new_ppar<<" pper " << pper << " lamda " <<lamda*Constants::R2D<< " alpha "<< alpha*Constants::R2D << " aeq " <<aeq*Constants::R2D<< " at time " << time ;
                single.escaped = true ;
            }
            break;
        }

        //Next step:
        new_values_RK4(lamda, ppar, pper, alpha, l1, l2, l3, l4, m1, m2, m3, m4, o1, o2, o3, o4, p1, p2, p3, p4);
        
        //Find aeq from alpha and lamda.
        aeq = asin(sin(alpha)*sqrt(Bmag_dipole(0)/Bmag_dipole(lamda))); 
        
        time  = time + Constants::h; 
        i++;  
       
		//To save any states:
		//single.save_state( p, lamda, alpha, aeq, ppar, pper, time);
        std::cout<<"\n\nalpha "<<alpha*Constants::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlamda " <<lamda*Constants::R2D<< "\naeq "<<aeq*Constants::R2D;

    }
    

    //Save last state to return values and continue the simulation with wave (if needed). 
    single.lamda_end = lamda;
    single.ppar_end  = ppar;
    single.pper_end  = pper;
    single.alpha_end = alpha;
    single.aeq_end   = aeq;
    single.time_end  = time;
    single.eta_end   = Constants::eta0;
    //single.zeta_end = zeta;
    //single.upar_end = upar;
    //single.uper_end = uper;
    //single.Ekin_end = Ekin;
}




