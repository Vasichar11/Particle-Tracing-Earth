#include "headers/time_rates.h"

//Equation 1: parallel(to B) speed uz=dz/dt, variation with time
real z_rk(real ppar_tmp, real gama)
{			
    real krk=ppar_tmp/(gama*Constants::m_e);
    return krk;
}

//Equation 2: parallel momentum dpz/dt, variation with time
real p_par_rk(real pper_tmp,real eta_tmp, real kz, real w_h, real dwh_ds, real gama, real wtau_sq)
{  
    real lrk;

    if(std::isnan(wtau_sq*Constants::m_e/kz))  //for no interaction to be valid.
    {
        lrk=-(1/(Constants::m_e*gama))*((pper_tmp*pper_tmp)/(2*w_h))*dwh_ds;  //((wtau_sq*m_e)/kz)*sin(eta_tmp)=0.
    }
    else
    {
        //Borntik thesis 2.24a
        lrk=((wtau_sq*Constants::m_e)/kz)*sin(eta_tmp)-(1/(Constants::m_e*gama))*((pper_tmp*pper_tmp)/(2*w_h))*dwh_ds;
    }    
    
    return lrk;
}

//Equation 3: differentiate perpendicular momentum, variation with time
real p_per_rk(real ppar_tmp,real pper_tmp,real eta_tmp, real w_h, real dwh_ds, real gama, real w1, real w2, real R1, real R2, real beta)
{
    //Borntik thesis 2.24b			
    real mrk=-(w1*((ppar_tmp/gama)+(Constants::m_e*R1))*jn((Constants::m_res)-1,beta)
			 - w2*((ppar_tmp/gama)-(Constants::m_e*R2))*jn((Constants::m_res)+1,beta))*sin(eta_tmp)+((1/(Constants::m_e*gama))*(pper_tmp*ppar_tmp)/(2*w_h))*dwh_ds;
    return mrk; 
}

//Equation 4: eta: angle between BwR and u_per, variation with time
real eta_rk(real kz,real ppar_tmp, real w_h, real gama)
{
    //Borntik thesis 2.24c
    real nrk=(((Constants::m_res)*w_h)/gama)-Constants::w_wave-kz*(ppar_tmp/(Constants::m_e*gama));
    return nrk;
}   

//Equation 5: // lamda variation with time
real lamda_rk(real ppar_tmp,real lamda_tmp,real gama)
{
    //Stelios eq (dz/dlamda=L_shell*Re*sqrt(1+3*sin(lamda)*sin(lamda))*cos(lamda))
    real ork=ppar_tmp/(gama*Constants::m_e*Constants::L_shell*Constants::Re*sqrt(1+3*sin(lamda_tmp)*sin(lamda_tmp))* cos(lamda_tmp));
    return ork;
}

//Equation 6:  P.A variation with time
real alpha_rk(real pper_tmp,real alpha_tmp,real eta_tmp,real kz, real w_h, real dwh_ds, real gama, real wtau_sq)
{
    real prk;

    if(std::isnan(Constants::m_e*wtau_sq/(kz*pper_tmp))) //for no interaction to be valid.
    {
        prk = ((1/(Constants::m_e*gama))*(pper_tmp/(2*w_h)))*dwh_ds;
    }
    else
    {   
        //Bortnik thesis 2.26
        prk=(-((Constants::m_e*wtau_sq/(kz*pper_tmp))*(1+((cos(alpha_tmp)*cos(alpha_tmp))/(Constants::m_res*(w_h/Constants::w_wave)-1))))*sin(eta_tmp)+
            ((1/(Constants::m_e*gama))*(pper_tmp/(2*w_h)))*dwh_ds);
    }

    return prk;
}
    
//Equation 7: Equatorial P.A variation with time.
real aeq_rk(real ppar_tmp, real pper_tmp, real alpha_tmp, real eta_tmp, real aeq_tmp, real kappa, real gama, real Bw_out)  
{
    real p_mag  = sqrt(ppar_tmp*ppar_tmp+pper_tmp*pper_tmp);
    real aeq_rk;

    if(kappa==0)    //(Constants::w_wave/kappa) goes to inf.
    {
        aeq_rk=0; //when no interaction
    }
    else
    {
        aeq_rk = ((Constants::q_e*Bw_out)/(pow(p_mag,2)))*(tan(aeq_tmp)/tan(alpha_tmp))*(((Constants::w_wave/kappa)-(ppar_tmp/(gama*Constants::m_e)))*ppar_tmp-(pow(pper_tmp,2)/(gama*Constants::m_e)))*sin(eta_tmp);
    }

    return aeq_rk;
}
