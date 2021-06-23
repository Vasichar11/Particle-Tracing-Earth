//Equation 1: $\frac{dz}{dt}$ #In[12]: //parallel(to B) speed uz=dz/dt, variation with time
real z_rk(real ppar_tmp, real m_e, real gama){			
    real krk=ppar_tmp/(gama*m_e);
return krk;
}

//Equation 2: $\frac{dp_z}{dt}$ #In[13]: //parallel momentum dpz/dt, variation with time
real p_par_rk(real pper_tmp,real eta_tmp, real m_e, real kz, real w_h, real dwh_ds, real gama, real wtau_sq){  
    //Borntik thesis 2.24a
    real lrk=((wtau_sq*m_e)/kz)*sin(eta_tmp)-(1/(m_e*gama))*((pper_tmp*pper_tmp)/(2*w_h))*dwh_ds;  
    
return lrk;
}

//Equation 3: $\frac{dp_{\perp}}{dt}$ #In[14]: //differentiate perpendicular momentum, variation with time
real p_per_rk(real ppar_tmp,real pper_tmp,real eta_tmp, real m_e, real w_h, real dwh_ds, real gama, real w1, real w2, real R1, real R2, real beta){
    //Borntik thesis 2.24b			
    real mrk=-(w1*((ppar_tmp/gama)+(m_e*R1))*jn((Constants::m_res)-1,beta) //replace m_res with 1
			 - w2*((ppar_tmp/gama)-(m_e*R2))*jn((Constants::m_res)+1,beta))*sin(eta_tmp)+((1/(m_e*gama))*(pper_tmp*ppar_tmp)/(2*w_h))*dwh_ds;
return mrk; //parentheseis?
}

//Equation 4: $ \frac{d\eta}{dt}$ #In[15]:	// eta: angle between BwR and u_per, variation with time
real eta_rk(real m_e,real w_wave,real kz,real ppar_tmp, real w_h, real gama){
    //Borntik thesis 2.24c
    real nrk=(((Constants::m_res)*w_h)/gama)-w_wave-kz*(ppar_tmp/(m_e*gama));
return nrk;
}   

//Equation 5: $ \frac{\partial \lambda}{\partial t}$ #In[16]:   // lamda variation with time
real lamda_rk(real ppar_tmp,real lamda_tmp,real m_e, real gama){
    //Stelios eq (dz/dlamda=L_shell*Re*sqrt(1+3*sin(lamda)*sin(lamda))*cos(lamda))
    real ork=ppar_tmp/(gama*m_e*Constants::L_shell*Constants::Re*sqrt(1+3*sin(lamda_tmp)*
                    sin(lamda_tmp))*cos(lamda_tmp));
return ork;
}

//Equation 6: $\frac{d\alpha}{dt}$ # In[17]:		// P.A variation with time
real alpha_rk(real pper_tmp,real alpha_tmp,real eta_tmp,real m_e, real kz, real w_h, real w, real dwh_ds, real gama, real wtau_sq){
    //Bortnik thesis 2.26

    real prk=(-((m_e*wtau_sq/(kz*pper_tmp))*(1+((cos(alpha_tmp)*cos(alpha_tmp))/(Constants::m_res*(w_h/w)-1))))*sin(eta_tmp)+
            ((1/(m_e*gama))*(pper_tmp/(2*w_h)))*dwh_ds);

return prk;
}