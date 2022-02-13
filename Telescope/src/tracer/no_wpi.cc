#include "headers/no_wpi.h"

//Adiabatic motion.
void no_wpi(int p, Particles &single, Telescope &ODPT, Particles &particle_state)
{
    real lamda    =  particle_state.lamda.front();
    real ppar     =  particle_state.ppar.front(); 
    real pper     =  particle_state.pper.front(); 
    real alpha    =  particle_state.alpha.front(); 
    real aeq      =  particle_state.aeq.front(); 
    real time     =  particle_state.time.front();
    //real zeta     =  particle_state.zeta.front(); 
    //real upar     =  particle_state.upar.front(); 
    //real uper     =  particle_state.uper.front();

    //Declare function's variables. Once for each particle. When parallel, declare xcore times?
    real new_lamda,new_ppar;
    real w_h, dwh_ds, Bmag, p_mag, gama;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,o1,o2,o3,o4,p1,p2,p3,p4;
    bool trapped=1;

    //Objects for each specie.
    Species electron(Constants::m_e,  Constants::q_e, 1); 
    Species oxygen  (Constants::m_O,  Constants::q_i, 0.006); 
    Species hydrogen(Constants::m_H,  Constants::q_i, 0.94); 
    Species helium  (Constants::m_He, Constants::q_i, 0.054);
    
	//std::cout.precision(64);			//Output 16 decimal precise
	//std::cout<<std::scientific;		//For e notation representation

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
            trapped = 0;
        }

        //Check Precipitation:
        new_ppar = ppar + (Constants::h/6)*(l1+2*l2+2*l3+l4);
        if(!trapped && (ppar*new_ppar<0) ) //Would bounce if ppar is about to change sign.
        {   
            //To save states of precipitating particles:
            #pragma omp critical //Only one processor should write at a time. Otherwise there is a chance of 2 processors writing in the same spot.
            {   
                single.save_state(p, lamda, alpha, aeq, ppar, pper, time);
                std::cout<<"\n\nParticle "<<p<<" escaped with ppar "<<ppar<< " new_ppar would be "<<new_ppar<<" pper " << pper << " lamda " <<lamda*Constants::R2D<< " alpha "<< alpha*Constants::R2D << " aeq " <<aeq*Constants::R2D<< " at time " << time ;
                
                //This particle won't continue for the WPI simulation.
                //Erase it's vectors. 
                particle_state.lamda.erase(particle_state.lamda.begin());
                particle_state.ppar.erase(particle_state.ppar.begin());
                particle_state.pper.erase(particle_state.pper.begin());
                particle_state.alpha.erase(particle_state.alpha.begin());
                particle_state.aeq.erase(particle_state.aeq.begin());
                particle_state.time.erase(particle_state.time.begin());
                particle_state.zeta.erase(particle_state.zeta.begin());
                particle_state.upar.erase(particle_state.upar.begin());
                particle_state.uper.erase(particle_state.uper.begin());
                particle_state.Ekin.erase(particle_state.Ekin.begin());
                particle_state.eta.erase(particle_state.eta.begin());
                particle_state.deta_dt.erase(particle_state.deta_dt.begin());
                particle_state.M_adiabatic.erase(particle_state.M_adiabatic.begin());
            }
            break;
        }

        //Next step:
        new_values_RK4(lamda, ppar, pper, alpha, l1, l2, l3, l4, m1, m2, m3, m4, o1, o2, o3, o4, p1, p2, p3, p4);
        time  = time + Constants::h; 
        i++;  

		//To save states:
		//single.save_state(aeq,alpha,lamda,deta_dt,time);
        //std::cout<<"\n\nalpha "<<alpha << "\nppar "<< ppar<< "\npper " << pper << "\nlamda " <<lamda<< "\naeq " <<aeq ;

        //Stop at equator:
        //if(eql_dstr[p].lamda.at(i)>0) {	
        //	break;}	
    }
    

    //Save last state to return values and continue the simulation with wave (if needed). 
    particle_state.lamda.front() = lamda;
    particle_state.ppar.front()  = ppar;
    particle_state.pper.front()  = pper;
    particle_state.alpha.front() = alpha;
    particle_state.aeq.front()   = aeq;
    particle_state.time.front()  = time;
    particle_state.eta.front()   = Constants::eta0;
    //particle_state.zeta.front() = zeta;
    //particle_state.upar.front() = upar;
    //particle_state.uper.front() = uper;
    //particle_state.Ekin.front() = Ekin;
}




