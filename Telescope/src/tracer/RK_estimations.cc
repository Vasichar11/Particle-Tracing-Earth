#include "headers/RK_estimations.h"


//Function to approximate all values of Runge Kutta's block.
void new_values_RK4(real &lamda, real &ppar, real &pper, real &eta, real &alpha, real &aeq, real l1, real l2, real l3, real l4, real m1, real m2, real m3, real m4, real n1, real n2, real n3, real n4, real o1, real o2, real o3, real o4, real p1, real p2, real p3, real p4, real q1, real q2, real q3, real q4)
{
    lamda  +=  (Constants::h/6)*(o1+2*o2+2*o3+o4);
    ppar   +=  (Constants::h/6)*(l1+2*l2+2*l3+l4);
    pper   +=  (Constants::h/6)*(m1+2*m2+2*m3+m4);
    eta    +=  (Constants::h/6)*(n1+2*n2+2*n3+n4);
    alpha  +=  (Constants::h/6)*(p1+2*p2+2*p3+p4);
    aeq    +=  (Constants::h/6)*(q1+2*q2+2*q3+q4);
    //zeta    =  zeta   +  (Constants::h/6)*(k1+2*k2+2*k3+k4);
    //deta_dt =            (Constants::h/6)*(n1+2*n2+2*n3+n4);
    //upar    =  ppar   /  (Constants::m_e*gama);
    //uper    =  pper   /  (Constants::m_e*gama);
    //p_mag = sqrt((ppar*ppar)+(pper*pper));
    //gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
    //Ekin = ((gama-1)*Constants::m_e*Constants::c*Constants::c)*6.2415e15; 
    //B_lam    =  Bmag_dipole(lamda);    
    //M_adiabatic = (pper*pper)/(2*Constants::m_e*B_lam); 
}

//Overloaded for noWPI 
//-> without eta,aeq estimation
void new_values_RK4(real &lamda, real &ppar, real &pper, real &alpha, real l1, real l2, real l3, real l4, real m1, real m2, real m3, real m4, real o1, real o2, real o3, real o4, real p1, real p2, real p3, real p4)
{
    lamda  +=  (Constants::h/6)*(o1+2*o2+2*o3+o4);
    ppar   +=  (Constants::h/6)*(l1+2*l2+2*l3+l4);
    pper   +=  (Constants::h/6)*(m1+2*m2+2*m3+m4);
    alpha  +=  (Constants::h/6)*(p1+2*p2+2*p3+p4);
    //zeta    =  zeta   +  (Constants::h/6)*(k1+2*k2+2*k3+k4);
    //deta_dt =            (Constants::h/6)*(n1+2*n2+2*n3+n4);
    //upar    =  ppar   /  (Constants::m_e*gama);
    //uper    =  pper   /  (Constants::m_e*gama);
    //p_mag = sqrt((ppar*ppar)+(pper*pper));
    //gama = sqrt((p_mag*p_mag*Constants::c*Constants::c)+(Constants::m_e*Constants::m_e*Constants::c*Constants::c*Constants::c*Constants::c))/(Constants::m_e*Constants::c*Constants::c);
    //Ekin = ((gama-1)*Constants::m_e*Constants::c*Constants::c)*6.2415e15; 
    //B_lam    =  Bmag_dipole(lamda);    
    //M_adiabatic = (pper*pper)/(2*Constants::m_e*B_lam); 
}