//Overloaded Functions for txt file creation (passing const by reference, performance should be ok )
void outfile(const std::vector<real> &vect1, const std::string& filename) 
{	
	std::ofstream output(filename);

	for(int i=0; i<vect1.size(); i++)
 	{
 		output << std::scientific << std::setprecision(16) << vect1.at(i)  << "\n" ; 
 	} 

}
void outfile(const std::vector<real> &vect1, const std::vector<real> &vect2, const std::string& filename) 
{	std::ofstream output(filename);
	for(int i=0; i<vect1.size(); i++)
 	{output << std::scientific << std::setprecision(16) << vect1.at(i) << " " << vect2.at(i) << "\n" ; } }
void outfile(const std::vector<real> &vect1, const std::vector<real> &vect2, const std::vector<real> &vect3, const std::string& filename) 
{	std::ofstream output(filename);
	for(int i=0; i<vect1.size(); i++)
 	{output << std::scientific << std::setprecision(16) << vect1.at(i) << " " << vect2.at(i)<< " " << vect3.at(i) << "\n" ; } }
void outfile(const std::vector<real> &vect1, const std::vector<real> &vect2, const std::vector<real> &vect3, const std::vector<real> &vect4, const std::string& filename)
{	std::ofstream output(filename);
	for(int i=0; i<vect1.size(); i++)
 	{output << std::scientific << std::setprecision(16) << vect1.at(i) << " " << vect2.at(i)<< " " << vect3.at(i) << " " << vect4.at(i) << "\n" ; } }
void outfile(const std::vector<real> &vect1, const std::vector<real> &vect2, const std::vector<real> &vect3, const std::vector<real> &vect4, const std::vector<real> &vect5, const std::string& filename) 
{	std::ofstream output(filename);
	for(int i=0; i<vect1.size(); i++)
 	{output << std::scientific << std::setprecision(16) << vect1.at(i) << " " << vect2.at(i)<< " " << vect3.at(i) << " " << vect4.at(i) << " " << vect5.at(i) << "\n" ; } }
void outfile(const std::vector<real> &vect1, const std::vector<real> &vect2, const std::vector<real> &vect3, const std::vector<real> &vect4, const std::vector<real> &vect5, const std::vector<real> &vect6, const std::string& filename) 
{	std::ofstream output(filename);
	for(int i=0; i<vect1.size(); i++)
 	{output << std::scientific << std::setprecision(16) << vect1.at(i) << " " << vect2.at(i)<< " " << vect3.at(i) << " " << vect4.at(i) << " " << vect5.at(i) << " " << vect6.at(i) << "\n" ; } }
void outfile(const std::vector<real> &vect1, const std::vector<real> &vect2, const std::vector<real> &vect3, const std::vector<real> &vect4, const std::vector<real> &vect5, const std::vector<real> &vect6, const std::vector<real> &vect7, const std::string& filename) 
{	std::ofstream output(filename);
	for(int i=0; i<vect1.size(); i++)
 	{output << std::scientific << std::setprecision(16) << vect1.at(i) << " " << vect2.at(i)<< " " << vect3.at(i) << " " << vect4.at(i) << " " << vect5.at(i) << " " << vect6.at(i) << " " << vect7.at(i) << "\n" ; } }



//Output object's variable members 
void outfile(const std::vector<Particles>& dstr, const std::string& filename)
{
	std::ofstream output(filename);
	
	for(int i=0; i<dstr[0].lamda.size(); i++)	// Line by line seperating states of time.
	{
		//For each particle in the population, variables are written between spaces.

		for(int p=0; p<dstr.size(); p++)									//eql_dstr[0].lamda(i) eql_dstr[1].lamda(i) ...
			output << std::setprecision(16) << dstr[p].lamda.at(i) << " " ;	//All lamdas of particles, space
		//for(int p=0; p<dstr.size(); p++)									
		//	output << std::setprecision(16) << dstr[p].zeta.at(i) << " " ;
		for(int p=0; p<dstr.size(); p++)									//... eql_dstr[0].aeq(i) eql_dstr[1].aeq(i) ...
			output << std::setprecision(16) << dstr[p].aeq.at(i) << " " ;	//All aeqs of particles, space...
		for(int p=0; p<dstr.size(); p++)	
			output << std::setprecision(16) << dstr[p].aeqsu.at(i) << " " ;
		for(int p=0; p<dstr.size(); p++)	
			output << std::setprecision(16) << dstr[p].alpha.at(i) << " " ;
		//for(int p=0; p<dstr.size(); p++)	
		//	output << std::setprecision(16) << dstr[p].eta.at(i) << " " ;
		//for(p=0; p<dstr.size(); p++)	
		//	output << std::setprecision(16) << dstr[p].upar.at(i) << " " ;
		//for(p=0; p<dstr.size(); p++)	
		//	output << std::setprecision(16) << dstr[p].uper.at(i) << " " ;
		//for(p=0; p<dstr.size(); p++)	
		//	output << std::setprecision(16) << dstr[p].ppar.at(i) << " " ;
		//for(p=0; p<dstr.size(); p++)	
		//	output << std::setprecision(16) << dstr[p].pper.at(i) << " " ;
		//for(p=0; p<dstr.size(); p++)	
		//	output << std::setprecision(16) << dstr[p].time.at(i) << " " ;	

		output << "\n" ;	//now next line        							 // i++
	}

}