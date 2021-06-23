//Temporal, there are already algorithms for csv reading...

//Function. Read filename. Return Vector of pairs <column name, column vector values>
std::vector<std::pair<std::string, std::vector<double>>> rcsv(std::string filename)
{	
	std::vector<std::pair<std::string, std::vector<double>>> output;	//Output vector
	std::ifstream file(filename);								  		//Input file stream, the csv

	if(file.is_open() && file.good() ) 				//Check if file is open/raising any errors
	{
		std::string line, column;
		std::getline(file, line); 					//Gets the first line of the file, i.e names of variables
		std::stringstream s_stream(line);			//As_streamociates line string with a stream we can read from
		while(std::getline(s_stream,column, ',')) 	//Extracts characters from s_stream, stores in column, until "," is met
		{	//push back pairs into the output
			output.push_back({column,std::vector<double> {}});
		}
//here, we yet have, an output pair with variable names and no values
		
		while(std::getline(file,line))			 	    	 //gets current line of the file
		{
			std::stringstream s_stream(line); 				 //Currect line stringstream
			int col = 0 ; 									 //Tracking columns
			double value;
			while(s_stream>>value)							 //Getting from the current stringstream the values
			{	
				output.at(col).second.push_back(value);		 //Pushing them back at current column
				if(s_stream.peek() ==',') s_stream.ignore(); //When comma is met, ignore it, continue. peek reads without extracting
				col++;
			}
		}
	}
	else throw std::runtime_error("Cant open");
	
	file.close();
	return output;
}
