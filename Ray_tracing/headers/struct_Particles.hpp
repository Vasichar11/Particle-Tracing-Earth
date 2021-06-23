//Struct for particle's state throughout the iterations
struct Particles
{ 				
	std::vector<real> lamda , zeta, uper , upar, ppar, pper, alpha, aeq, aeqsu, eta, time; 
	//Constructor for initial state, called once when object is created.
	Particles(real first_lamda, real first_zeta, real first_upar, real first_uper, real first_ppar, real first_pper, real first_alpha, real first_aeq, real first_aeqsu, real first_time)
	{ 
		lamda.push_back(first_lamda);      		//Make alternative constructor.			
		zeta.push_back(first_zeta);
		upar.push_back(first_upar);
		uper.push_back(first_uper);	     
		ppar.push_back(first_ppar);		 
		pper.push_back(first_pper);
		alpha.push_back(first_alpha);
		aeq.push_back(first_aeq);
		aeqsu.push_back(first_aeq);
		//eta.push_back(first_eta);			
		time.push_back(first_time);
		//L_shell of particle will be added
	}   

	//Member function to push back initial gyrophases.
	void set_eta(real eta_distribution)
	{
		eta.push_back(eta_distribution);
	}

	//Member function to push_back new state
	void update_state(real new_lamda, real new_zeta, real new_uper, real new_upar, real new_ppar, real new_pper, real new_alpha, real new_aeq, real new_aeqsu, real new_eta, real new_time)
	{
		lamda.push_back(new_lamda);      				
		zeta.push_back(new_zeta);
		upar.push_back(new_upar);
		uper.push_back(new_uper);	     
		ppar.push_back(new_ppar);		 
		pper.push_back(new_pper);
		alpha.push_back(new_alpha);
		aeq.push_back(new_aeq);
		aeqsu.push_back(new_aeqsu);
		eta.push_back(new_eta);			
		time.push_back(new_time);
	}
	

};  	
