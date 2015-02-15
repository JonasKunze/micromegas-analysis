/*
	 data variables:
	 1D vector event->apv_id: 		ID of APV which shows hit (needed to select if hit is in X or Y)
	 1D vector event->mm_strip:		number of the strip which shows signal
	 1D vector event->apv_qmax: 	maximum charge of the strip
	 1D vector event->apv_tbqmax: 	time section of maximum charge (time = #TimeSlice * 25)

	 access element of 1D vector: eg. event->apv_id->at(i) gives strip data at position i in all vectors correspond to

	 2-D vector event->apv_q: 	matrix with full time characteristics of the signal

	 access element of 2D vector: event->apv_q->at(i).at(j) gives charge at time step j of strip event->apv_id->at(i)


	 useful code:
	 event->apv_q->size() gives size of this vector (size is the same for all 1D vector and first dimension of 2D vector)
	 event->apv_q->at(i).size() gives number of recorded timesteps

	 list:
	 type listName[numberofEntries];		//initialize list, number of entries cannot be changed afterwards, type can be int, float, ...
	 float bla[8];
	 listName[i];				//access element i of the list

	 vector:
	 vector<type> vectorName; 		//initialise vector called vectorName, compared to list no fixed lenght, type can be int, float, ...
	 vectorName.push_back(value);		//append value to vector
	 vectorName.size();			//returns number of elements in the vector
	 vectorName.at(i);			//returns element i

	 if:
	 if(fun>0){
	 cout << "you should work" << endl;
	 }else{
	 cout << "take a break" << endl;
	 }

	 combined conditions:
	 condition1 && condition2: 	condition1 AND condition2
	 condition1 || condition2: 	condition1 OR condition2


	 loop:
	 for(int i=0; i<10; i++){
	 cout << "Micromegas are great" << endl;
	 }

	 initialize a histogram with eventnumber in the name:
	 stringstream nameOfStringstream;
	 nameOfStringstream.str("");
	 nameOfStringstream << eventNumber << "nameOfHistogram";  //eg: for eventNumber=10 name will be 10nameOfHistogram
	 general_mapHist2DEvent[nameOfStringstream.str()] = new TH2F(nameOfStringstream.str().c_str(),";label x-axis; label y-axis",number of bins in x,smallest value x ,largest value x, number of bins in y, smallest value y, largest value x );

	 Gaussian Fit to histogram:
	 1. initialise histogram and function to store the fit
	 TH1F *histName = new TH1F("histname",";label label x-axis; label y-axis", number of bins, smallest value, largest value);
	 TF1 *fitName;
	 2. fill histogram
	 histName->SetBinContent(number of bin, value);
	 3. fit gauss to histogram
	 histName->Fit("gaus","q");
	 fitName = histName->GetFunction("gaus");
	 4. access parameters of the fit
	 fitName->GetParameter(parameterNumber); 	//parameterNumber for Gaussion fit: 0 = amplitude, 1 = mean, 2 = standard deviation
	 fitName->GetParError(paramterNumber);
	 fitName->GetChisquare();
	 fitName->GetNDF();				//has to be cast to double to fill in tree: (double) fitName->GetNDF();

	 */
	/*TODO:
	 1. Create 2D Eventdisplay (strips:time:charge)
	 2. Find maximum charge and store value and corresponding strip, time step, etc.
	 3. Develop cuts to select good events (check with Eventdisplay)
	 4. Gaussian fits to charge distribution over strips at timestep with maximum charge
	 5. store values in designated histograms and trees
	 */