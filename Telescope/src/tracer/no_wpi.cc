#include "headers/no_wpi.h"

//Adiabatic motion.
void no_wpi(int p, Particles &single, Telescope &ODPT)
{
	//std::cout.precision(64);			//Output 16 decimal precise
	//std::cout<<std::scientific;		//For e notation representation

    //Assign first particle states.
    real lamda     =  single.lamda0;
    real ppar      =  single.ppar0;
    real pper      =  single.pper0;
    real alpha     =  single.alpha0;
    real aeq       =  single.aeq0;
    real time      =  single.time0;
    real Ekin      =  single.Ekin0;
    //real zeta      =  single.zeta0;
    //real upar      =  single.upar0;
    //real uper      =  single.uper0;
    //std::cout<<"\n\ntime " << time << "\nalpha "<<alpha*Constants::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlamda " <<lamda*Constants::R2D<< "\naeq "<<aeq*Constants::R2D;

    //Declare function's variables. Once for each particle. When parallel, declare xcore times?
    real new_lamda, new_ppar, new_aeq, new_Ekin;
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



        //Runge kutta 4 first estimations:
        //First make the next step increments(lamda,aeq,ppar) that are needed to characterize crossing-validity-trapping. 
        new_lamda = lamda + (Constants::h/6)*(o1+2*o2+2*o3+o4);
        //Adiabatic motion. Aeq stays the same or changes between two values(?)
        //Is this ok?
        int k;
        if(ppar<0) k=1;
        else k=0;
        new_aeq = pow(-1,k) * asin(sin(alpha)*sqrt(Bmag_dipole(0)/Bmag_dipole(new_lamda))) + k*M_PI; 
        //We should have new_aeq   = aeq; and new_Ekin  = Ekin; for adiabatic motion. In practice, aeq changes between two values. One for northward motion and one for southward while on equator. (θ and 180-θ)
        new_Ekin   = ((gama-1)*Constants::m_e*Constants::c*Constants::c)*6.2415e15;
        new_ppar  = ppar  + (Constants::h/6)*(l1+2*l2+2*l3+l4);


        //Check Validity:
        if(std::isnan(new_lamda*new_aeq*new_ppar))
        {
            single.nan = true;
            single.nan_state( p ); //Save the initial state (after noWPI) for the particle and the time that P.A turned out negative.
            std::cout<<"\nParticleV "<<p<<" nan";
            break; 
        }
        //Check Negative P.A:
        if(alpha<0 || aeq<0)
        {
            single.negative = true;
            single.negative_state( p ); //Save the initial state (after noWPI) for the particle and the time that P.A turned out negative.
            std::cout<<"\nParticleN "<<p<<" negative p.a";
            break;
        }
        //Check higher than 180 P.A:
        if(alpha>M_PI)
        {
            single.high = true;
            single.high_state( p ); //Save the initial state (after noWPI) for the particle and the time that P.A turned out higher than 180.
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
        
        
        //Runge kutta 4 estimations:
        lamda   = new_lamda;
        aeq     = new_aeq;
        ppar    = new_ppar;
        Ekin    = new_Ekin;
        //Rest increments for the next step:
        alpha  +=  (Constants::h/6)*(p1+2*p2+2*p3+p4);
        pper   +=  (Constants::h/6)*(m1+2*m2+2*m3+m4);
        
        time  = time + Constants::h; 
        i++;  
       
        //std::cout<<"\n\ntime " << time <<" latitude " << lamda*Constants::R2D<< " Ekin " << Ekin << "\nalpha "<<alpha*Constants::R2D << " aeq0 " << single.aeq0 * Constants::R2D << "\nppar "<< ppar<< "\npper " << pper << "\nlamda " <<lamda*Constants::R2D<< "\naeq "<<aeq*Constants::R2D ;

    }
    

    //Save last state to return values and continue the simulation with wave (if needed). 
    single.lamda00 = lamda;
    single.ppar00  = ppar;
    single.pper00  = pper;
    single.alpha00 = alpha;
    single.aeq00   = aeq;
    single.time00  = time;
    single.Ekin00  = Ekin;
    single.eta00   = Constants::eta0; //Particle "gyrophase" initialization 
    //single.zeta00 = zeta;
    //single.upar00 = upar;
    //single.uper00 = uper;
}




