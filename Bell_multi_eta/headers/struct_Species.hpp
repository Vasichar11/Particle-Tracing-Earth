struct Species{     	  	      //Struct for mass & charge data of species				  
	real mass;		   	  	 	  //    and for computing densities & wps,wc
	real charge;
	real n_factor;        	 	  //For now, n_factor: e->1 , H->0.94, He->0.0054, O->0.006
	real density(real l)		  //Returns specie's lamda dependent density(estimation)
	{
		real clat=cos(l); 				
		real ns=(Constants::ne_0)*pow(clat,-4)*n_factor;  
		return ns;
	}
	real wps(real ns){ return ns*pow(charge,2)/(mass*(Constants::eps0)); } //Returns specie's plasma frequency squared
	real wc(real B0mag){ return (charge * B0mag)/mass; } 				   //Returns specie's cyclotron frequency
};