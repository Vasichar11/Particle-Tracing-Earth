#include "headers/adiabatic_motion.h"

//Adiabatic motion.
void adiabatic_motion(int64_t track_pop,int p, Particles &single, Telescope &ODPT)
{
	std::cout<<"\rParticle "<<p<<" is bouncing"<<std::flush;
    
    
    real lamda =  single.lamda.at(0);
    real zeta  =  single.zeta.at(0); 
    real ppar  =  single.ppar.at(0); 
    real pper  =  single.pper.at(0); 
    real eta   =  single.eta.at(0); 
    real alpha =  single.alpha.at(0); 
    real aeq   =  single.aeq.at(0); 
    real upar  =  single.upar.at(0); 
    real uper  =  single.uper.at(0); 
    real time  =  single.time.at(0);

    //Declare function's variables. Once for each particle. When parallel, declare xcore times?
    real new_lamda;
    real wc_e,w_h, dwh_ds, Bmag, p_mag, gama;
    real k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4,n1,n2,n3,n4,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4;

    //Objects for each specie. Used inside "motion" function. //Change that, it's called once for each particle...
    Species electron(Constants::m_e,  Constants::q_e, 1); 
    Species oxygen  (Constants::m_O,  Constants::q_i, 0.006); 
    Species hydrogen(Constants::m_H,  Constants::q_i, 0.94); 
    Species helium  (Constants::m_He, Constants::q_i, 0.054);
    std::vector <std::vector<real>> time_sim  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );

    //Make vectors to save data throughout the Runge Kutta only if needed.(0) 
    //Create vectors instead of big arrays to allocate in heap.
    //Each row(population) has a vector which has cols(Nsteps+1) number of elements, each element initialized to 0.

    //std::vector <std::vector<real>> Ekin      (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) ); 
    //std::vector <std::vector<real>> wh_out    (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> gama_out  (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>> deta_dt   (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>>dwh_dt_out (track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
    //std::vector <std::vector<real>>B_earth_out(track_pop, std::vector<real> (Constants::Nsteps + 1, 0) );
    
	std::cout.precision(8);			//Output 16 decimal precise
	std::cout<<std::scientific;		//For e notation representation

    int i=0;
    
    while(i<Constants::Nsteps) 
    {
        
        Bmag=Bmag_dipole(lamda);
        //Cyclotron frequencies.
        wc_e = w_h = electron.wc(Bmag);
        //dwh_dt.
        dwh_ds=dwh_dsf(w_h,lamda);
        //From bell parameters function(temp for no interaction).
        p_mag = sqrt(ppar*ppar+pper*pper);
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //Check print parameters
        //std::cout<<"\n" << "\nwc_O " <<wc_O << "\ngama " << gama;     
        //Keep these values in arrays only if needed.(1) 
        //B_earth_out[p][i]=Bmag;  //wh_out[p][i] = w_h; //dwh_dt_out[p][i] = dwh_ds; //gama_out[p][i]=gama;
        //RK step-1//#################################################################################################################################################################################################################################################################################################
        k1=z_rk(ppar,gama);
        l1=p_par_rk(pper,eta,w_h,dwh_ds,gama);
        m1=p_per_rk(ppar,pper,eta,w_h,dwh_ds,gama);
        n1=eta_rk(ppar,w_h,gama);
        o1=lamda_rk(ppar,lamda,gama);
        p1=alpha_rk(pper,w_h,dwh_ds,gama);
        q1=0; 
        //std::cout<<"\n" << "k1 " << k1 << "\nl1 " <<l1 << "\nm1 " << m1 << "\nn " << n1<< "\no1 " << o1 << "\np1 " << p1 << "\nq1 " << q1 <<"\n";	
        
        
        Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o1);
        wc_e = w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o1);
        p_mag = sqrt((ppar+0.5*(Constants::h)*l1)*(ppar+0.5*(Constants::h)*l1)+(pper+0.5*(Constants::h)*m1)*(pper+0.5*(Constants::h)*m1));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-2//#################################################################################################################################################################################################################################################################################################
        k2=z_rk(ppar+0.5*(Constants::h)*l1,gama);  
        l2=p_par_rk(pper+0.5*(Constants::h)*m1,eta+0.5*(Constants::h)*n1,w_h,dwh_ds,gama);
        m2=p_per_rk(ppar+0.5*(Constants::h)*l1,pper+0.5*(Constants::h)*m1,eta+0.5*(Constants::h)*n1,w_h,dwh_ds,gama);
        n2=eta_rk(ppar+0.5*(Constants::h)*l1,w_h,gama);
        o2=lamda_rk(ppar+0.5*(Constants::h)*l1,lamda+0.5*(Constants::h)*o1,gama);
        p2=alpha_rk(pper+0.5*(Constants::h)*m1,w_h,dwh_ds,gama);
        q2=0;
        //std::cout<<"\n" << "k2 " << k2 << "\nl2 " <<l2 << "\nm2 " << m2 << "\nn2 " << n2<< "\no2 " << o2 << "\np2 " << p2 << "\n";
        

        Bmag=Bmag_dipole(lamda+0.5*(Constants::h)*o2);
        wc_e = w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda+0.5*(Constants::h)*o2);
        p_mag = sqrt((ppar+0.5*(Constants::h)*l2)*(ppar+0.5*(Constants::h)*l2)+(pper+0.5*(Constants::h)*m2)*(pper+0.5*(Constants::h)*m2));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-3//#################################################################################################################################################################################################################################################################################################
        k3=z_rk(ppar+0.5*(Constants::h)*l2,gama);  
        l3=p_par_rk(pper+0.5*(Constants::h)*m2,eta+0.5*(Constants::h)*n2,w_h,dwh_ds,gama);
        m3=p_per_rk(ppar+0.5*(Constants::h)*l2,pper+0.5*(Constants::h)*m2,eta+0.5*(Constants::h)*n2,w_h,dwh_ds,gama);
        n3=eta_rk(ppar+0.5*(Constants::h)*l2,w_h,gama);
        o3=lamda_rk(ppar+0.5*(Constants::h)*l2,lamda+0.5*(Constants::h)*o2,gama);
        p3=alpha_rk(pper+0.5*(Constants::h)*m2,w_h,dwh_ds,gama);
        q3=0;
        //std::cout<<"\n" << "k3 " << k3 << "\nl3 " <<l3 << "\nm3 " << m3 << "\nn3 " << n3<< "\no3 " << o3 << "\np3 " << p3 << "\n";


        Bmag=Bmag_dipole(lamda+(Constants::h)*o3);
        wc_e = w_h = electron.wc(Bmag);
        dwh_ds=dwh_dsf(w_h,lamda+(Constants::h)*o3);
        p_mag = sqrt((ppar+(Constants::h)*l3)*(ppar+(Constants::h)*l3)+(pper+(Constants::h)*m3)*(pper+(Constants::h)*m3));
        gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
        //RK step-4//#################################################################################################################################################################################################################################################################################################																								
        k4=z_rk(ppar+(Constants::h)*l3,gama);  
        l4=p_par_rk(pper+(Constants::h)*m3,eta+(Constants::h)*n3,w_h,dwh_ds,gama);
        m4=p_per_rk(ppar+(Constants::h)*l3,pper+(Constants::h)*m3,eta+(Constants::h)*n3,w_h,dwh_ds,gama);
        n4=eta_rk(ppar+(Constants::h)*l3,w_h,gama);
        o4=lamda_rk(ppar+(Constants::h)*l3,lamda+(Constants::h)*o3,gama);
        p4=alpha_rk(pper+(Constants::h)*m3,w_h,dwh_ds,gama);
        q4=0; 
        //std::cout<<"\n" << "k4 " << k4 << "\nl4 " <<l4 << "\nm4 " << m4 << "\nn " << n4<< "\no4 " << o4 << "\np4 " << p4 << "\n";

        

        //Approximate new lamda
        new_lamda = lamda + ((Constants::h)/6)*(o1+2*o2+2*o3+o4);
        
        //Check if NAN to break this particle. Why NAN ?
        if(std::isnan(new_lamda))
        {
            std::cout<<"\n\nParticle "<<p<<" with aeq0="<<" breaks.";
            //std::cout<<"\n"<< alpha << " " << zeta << " " << ppar<< " " << pper<< " " << eta << " " <<lamda<< " " <<aeq ;
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

        //deta_dt[p][i+1] = (Constants::h/6)*(n1+2*n2+2*n3+n4);

        //B_lam    =  Bmag_dipole(lamda);    
        //M_adiabatic = (pper*pper)/(2*Constants::m_e*B_lam); 

        //Go to next timestep
        time  = time + Constants::h; 
        
		//To save states:
		single.save_state(aeq,alpha,lamda,time);

        i++;  

        //Stop at equator
        //if(eql_dstr[p].lamda.at(i)>0) {	
        //	break;}	
    }
    //std::cout<<"\n\nzeta "<< zeta << "\nppar "<< ppar<< "\npper " << pper<< "\neta " << eta << "\nlamda " <<lamda<< "\nalpha "<< alpha << "\naeq " <<aeq ;

}




