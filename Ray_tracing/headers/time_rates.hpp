void slopes(real &k, real &l, real &m, real &n, real &o, real &p, const real ppar_tmp, const real pper_tmp, const real alpha_tmp, const real lamda_tmp, const real eta_tmp, const real Fpar, const real Fper, const real Ftheta, const real gama, const real w_h, const real dwh_ds, const real kz)
{	

	using namespace Constants;


	//Speed Parallel(to the static magnetic field Bo) time rate. 
	k = ppar_tmp/(gama*m_e);
	
	//Parallel momentum variation with time. Shouldn't it be f(time,par). Where's parallel momentum dependence?if not->Simpson's rule     
	l = Fpar*sin(eta_tmp)-(1/(m_e*gama))*((pper_tmp*pper_tmp)/(2*w_h))*dwh_ds;
	
	//Perpendicular momentum variation with time       
	m = Fper*sin(eta_tmp) + ((1/(m_e*gama))*(pper_tmp*ppar_tmp)/(2*w_h))*dwh_ds;
	
	//eta  variation with time #In[15]:     //eta: Angle between BwR and u_per, w_h local electron gyrofreq.
	n = ((m_res*Ftheta)/pper_tmp)*cos(eta_tmp)  +  ((m_res*w_h)/gama)  -  w_wave  -  kz*(ppar_tmp/(m_e*gama));
	
	//Lamda time rate.    
	o = ppar_tmp / (gama * m_e * L_shell * Re * sqrt(1+3*sin(lamda_tmp)*sin(lamda_tmp)) * cos(lamda_tmp));
	
	//P.A variation time rate.                   
	p = (-(Fpar/pper_tmp)*(1+(cos(alpha_tmp)*cos(alpha_tmp))/(Constants::m_res*(w_h/(gama*Constants::w_wave))-1))*sin(eta_tmp) + ((1/(m_e*gama))*(pper_tmp/(2*w_h)))*dwh_ds);



}