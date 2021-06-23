struct Waves{					 //Struct for wave field components
	real By;					 //Known magnetic field
	real w;						 //Angular frequency of the wave
	real f;						 //Frequency of the wave
	real th; 					 //Angle
	real th_d;
	std::vector<real> Bw,Bxw,Byw,Bzw,Ew,Exw,Eyw,Ezw;
	Waves(int64_t size) 
	{
		Bw.resize(size);
		Bxw.resize(size);
		Byw.resize(size);
		Bzw.resize(size);
		Ew.resize(size);
		Exw.resize(size);
		Eyw.resize(size);
		Ezw.resize(size);
	}
	void whistlers(real i,real mu,real P,real D,real S)
	{	//Compute field components
		real mu_sq = pow(mu,2); //mu*mu inaccurate
		real fac1 = P-(mu_sq*sin(th)*sin(th)) ;
		Bxw.at(i)=std::abs((-(D*fac1)/(P*(S-mu_sq)))*By); 
		Byw.at(i)=By;
		Bzw.at(i)=std::abs(((D*sin(th)*fac1)/(P*cos(th)*(S-(mu_sq)))*By));
		Exw.at(i)=std::abs((((Constants::c)*fac1)/(mu*P*cos(th))*By));
		Eyw.at(i)=std::abs(((D*(Constants::c)*fac1)/(mu*P*cos(th)*(mu_sq-S)))*By);
		Ezw.at(i)=std::abs((-((Constants::c)*mu*sin(th))/P)*By);
	}
};