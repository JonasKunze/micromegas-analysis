#include "CCommonIncludes.h"
#include "MMQuickEvent.h"
#include "MapFile.h"
#include "TGraphErrors.h"

using namespace std;

// Global Variables
MMQuickEvent *m_event;
map<string, TTree*> general_mapTree; 		//TTress
map<string, TH1F*> general_mapHist1D; 	//1D histogram of analysis for each run
map<string, TH2F*> general_mapHist2D;	//2D histogram of analysis for each run
map<string, TH2F*> general_mapHist2DEvent; 	//plot of 2D EventDisplay
map<string, TH1F*> general_mapPlotFit;		//plot of fits
map<string, TH2F*> general_mapCombined;		//combined Plots

Double_t m_TotalEventNumber;
vector<double> eventTimes;

//structure for trees
struct gauss_t {
	Double_t gaussXmean;
	Double_t gaussXmeanError;
	Double_t gaussXsigma;
	Double_t gaussXcharge;
	Double_t gaussXchi;
	Double_t gaussXdof;
	Double_t gaussXchiRed;
	Double_t gaussYmean;
	Double_t gaussYmeanError;
	Double_t gaussYsigma;
	Double_t gaussYcharge;
	Double_t gaussYchi;
	Double_t gaussYdof;
	Double_t gaussYchiRed;
	Int_t number;
};

struct maxi_t {
	Int_t maxXmean;
	Int_t maxXcharge;
	Int_t maxXcluster;
	Int_t maxYmean;
	Int_t maxYcharge;
	Int_t maxYcluster;
	Int_t number;
};

gauss_t gauss;
maxi_t maxi;

//set output path and name of output files
const string inPath = "/localscratch/praktikum/data/";		//Path of the Input
const string outPath = "/localscratch/praktikum/output/"; // Path of the Output
const string appendName = "";					// Name of single measurements
const string combinedPlotsFile = "";// Name of the file for the combined results of all runs

// Mapping of APV-Chips
const int APVIDMM_X0 = 5;
const int APVIDMM_X1 = 4;
const int APVIDMM_X2 = 6;
const int APVIDMM_Y0 = 0;
const int APVIDMM_Y1 = 1;
const int APVIDMM_Y2 = 2;

//number of strips in x and y
const int xStrips = 360;
const int yStrips = 360;

//Voltage range, needed for initialization of combined histograms
const int driftStart = 100;
const int driftEnd = 400;
const int driftSteps = 50;
const int ampStart = 500;
const int ampEnd = 550;
const int ampSteps = 25;

// functions to select if hit is in X or Y according to APV ID and mapping while data aquisition
bool isX(int id) {
	if (id == APVIDMM_X0 || id == APVIDMM_X1 || id == APVIDMM_X2) {
		return true;
	} else {
		return false;
	}
}

bool isY(int id) {
	if (id == APVIDMM_Y0 || id == APVIDMM_Y1 || id == APVIDMM_Y2) {
		return true;
	} else {
		return false;
	}
}

/**
 * Generates a new 3D-Histogram showing all measured charges in all strips of all times slices.
 * The Histogram will be stored at general_mapHist2DEvent[eventNumber+"nameOfHistogram"]
 */
TH2F* generateEventDisplay(MMQuickEvent* event, int eventNumber) {

	vector<vector<short> > chargeByStripByTime = *event->apv_q;
	unsigned int numberOfStrips = chargeByStripByTime.size();
	unsigned int numberOfTimeSlices = chargeByStripByTime[0].size();

	/*
	 * Generate a new 3D histogram for the event display (x=strip, y=timeslice, z=charge)
	 *
	 */

	// Generate the title of the histogram
	stringstream histoName;
	histoName.str("");
	histoName << eventNumber << "-Eventdisplay";

	/*
	 * Initialize a new root TH2F histogram with the right title, labels and bining:
	 *
	 * We'll have numberOfStrips x-bins going from 0 to numberOfStrips -1 and
	 * numberOfTimeSclies x-bins going from 0 to numberOfTimeSlices-1
	 */
	TH2F* eventDisplay = new TH2F(histoName.str().c_str(),
			";Strip Number; Time [25 ns]", numberOfStrips, 0,
			numberOfStrips - 1, numberOfTimeSlices, 0, numberOfTimeSlices - 1);

	// Store the new histogram in the global map
	general_mapHist2DEvent[histoName.str()] = eventDisplay;

	/*
	 * Fill eventDisplay with all measured charges at all strips for all times
	 */
	for (unsigned int stripNum = 0; stripNum != numberOfStrips; stripNum++) {
		for (unsigned int timeSlice = 0; timeSlice != numberOfTimeSlices;
				timeSlice++) {
			short charge = chargeByStripByTime[stripNum][timeSlice];
			eventDisplay->SetBinContent(stripNum/*x*/, timeSlice/*y*/,
					charge/*z*/);
		}
	}
	return eventDisplay;
}

void fitGauss(vector<short> chargeByStripAtMaxChargeTime, int eventNumber) {
	// Generate the title of the histogram
	stringstream histoName;
	histoName.str("");
	histoName << eventNumber << "-maxChargeDistribution";

	TH1F *maxChargeDistribution = new TH1F(histoName.str(), "; strip; charge",
			chargeByStripAtMaxChargeTime.size(), 0,
			chargeByStripAtMaxChargeTime.size() - 1);

	general_mapPlotFit[histoName.str()] = maxChargeDistribution;

	/*
	 * Fill the histogram
	 */
	for (unsigned int strip = 0; strip != chargeByStripAtMaxChargeTime.size();
			strip++) {
		maxChargeDistribution->SetBinContent(strip,
				chargeByStripAtMaxChargeTime[strip]);
	}

	maxChargeDistribution->Fit("gaus", "q");
	TF1* gaussFit = maxChargeDistribution->GetFunction("gaus");

	double mean = gaussFit->GetParameter(1); //parameterNumber for Gaussion fit: 0 = amplitude, 1 = mean, 2 = standard deviation
	double dMean = gaussFit->GetParError(1);

	double sigma = gaussFit->GetParameter(2); //parameterNumber for Gaussion fit: 0 = amplitude, 1 = mean, 2 = standard deviation
	double dSigma = gaussFit->GetParError(2);

//	gaussFit->GetChisquare();
//	gaussFit->GetNDF();
	double reducedChiSquare = gaussFit->GetChisquare() / gaussFit->GetNDF();
}

// analysis of single event: characteristics of event and gaussian fit
bool analyseMMEvent(MMQuickEvent *event, int eventNumber, int TRGBURST) {

	/*
	 data variables:
	 1D vector event->apv_id: 		ID of APV which shows hit (needed to select if hit is in X or Y)
	 1D vector event->mm_strip:		number of the strip which shows signal
	 1D vector event->apv_qmax: 	maximum charge of the strip
	 1D vector event->apv_tbqmax: 	time slice of maximum charge (time = timeSlice * 25)

	 access element of 1D vector: eg. event->apv_id->at(i) gives strip data at position i in all vectors correspond to

	 2-D vector event->apv_q: 	matrix with full time characteristics of the signal

	 access element of 2D vector: event->apv_q->at(i).at(j) gives charge at time step j of strip event->apv_id->at(i)


	 usefull code:
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

	//MMQuickEvent *event, int eventNumber, int TRGBURST
	vector<unsigned int> xOrY = *event->apv_id;
	vector<unsigned int> stripNumShowingSignal = *event->mm_strip;
	vector<short> maxChargeOfStrip = *event->apv_qmax;
	vector<short> timeSliceOfMaxChargeOfStrip = *event->apv_tbqmax;
	vector<vector<short> > chargeByStripByTime = *event->apv_q;

	/*
	 * 1. Created Eventdisplay
	 */
	generateEventDisplay(event, eventNumber);

	/*
	 * 2. Find maximum charge
	 */
	short maxCharge = maxChargeOfStrip[0];
	int stripWithMaxCharge = 0;
	int timeSliceOfMaxCharge = 0;

	for (unsigned int strip = 0; strip != maxChargeOfStrip.size(); strip++) {
		if (maxChargeOfStrip[strip] > maxCharge) {
			maxCharge = maxChargeOfStrip[strip];
			stripWithMaxCharge = strip;
			timeSliceOfMaxCharge = timeSliceOfMaxChargeOfStrip[strip];
		}
	}

	/*
	 * 4. Gaussian fits to charge distribution over strips at timestep with maximum charge
	 */
	vector<short> chargeByStripAtMaxChargeTime;
	for (unsigned int strip = 0; strip != chargeByStripByTime.size(); strip++) {
		vector<short> chargeByTime = chargeByStripByTime[strip];
		chargeByStripAtMaxChargeTime[strip] =
				chargeByTime[timeSliceOfMaxCharge];
	}

	fitGauss(chargeByStripAtMaxChargeTime, eventNumber);

//storage after procession
//Fill trees	(replace 1)
	if (/*condition to store the fit*/false) {
		gauss.gaussXmean = 1;
		gauss.gaussXmeanError = 1;
		gauss.gaussXsigma = 1;
		gauss.gaussXcharge = 1;
		gauss.gaussXchi = 1;
		gauss.gaussXdof = 1;
		gauss.gaussXchiRed = 1;
		gauss.gaussYmean = 1;
		gauss.gaussYmeanError = 1;
		gauss.gaussYsigma = 1;
		gauss.gaussYcharge = 1;
		gauss.gaussYchi = 1;
		gauss.gaussYdof = 1;
		gauss.gaussYchiRed = 1;
		gauss.number = eventNumber;

		maxi.maxXmean = 1;
		maxi.maxYmean = 1;
		maxi.maxXcharge = 1;
		maxi.maxYcharge = 1;
		maxi.maxXcluster = 1;
		maxi.maxYcluster = 1;
		maxi.number = eventNumber;

		general_mapTree["fits"]->Fill();

		general_mapHist2D["mmhitmapMaxCut"]->Fill(/*x value*/1,/*y value*/1);
		general_mapHist2D["mmhitmapGausCut"]->Fill(/*x value*/1,/*y value*/1);

	}

	if (/*condition to fill general histograms*/false) {
		if (/*condition to fill hitmap*/false) { // coincidence check between x and y signal (within 25ns)
			general_mapHist2D["mmhitmap"]->Fill(
			/*strip with maximum charge in X*/1,/*strip with maximum charge in Y*/
			1);
		}
		general_mapHist1D["mmchargex"]->Fill(/*maximum charge x*/1);
		general_mapHist1D["mmchargey"]->Fill(/*maximum charge y*/1);
		general_mapHist1D["mmhitx"]->Fill(/*strip x with maximum charge*/1);
		general_mapHist1D["mmhity"]->Fill(/*strip y with maximum charge*/1);
		general_mapHist1D["mmclusterx"]->Fill(
		/*number of x strips hit by one event*/1);
		general_mapHist1D["mmclustery"]->Fill(
		/*number of y strips hit by one event*/1);
		general_mapHist1D["mmtimex"]->Fill(/*time of maximum charge*/1);
		general_mapHist1D["mmtimey"]->Fill(/*time of maximum charege*/1);
		eventTimes.push_back(
				(double) event->time_s + (double) event->time_us / 1e6);

	}

	return true;
}

// Main Program
int main(int argc, char *argv[]) {

// map files to read different run of data in a row
// get data file name from MapFile.h
	MapFile* MicroMegas = new MapFile(inPath, outPath, appendName);
	map<string, TFile*> mapFile = MicroMegas->getFile();

//TRGBURST gives number of recorded timesteps (variable from data aquisition)
//timesteps = (TRGBURST+1)*3
	const int TRGBURST = 8;

//initialize file and histograms for combined output of all runs
	int numberOfXBins = (driftEnd - driftStart) / driftSteps + 1;
	double firstXBinValue = driftStart - 0.5 * driftSteps;
	double lastXBinValue = driftEnd + 0.5 * driftSteps;

	TFile* fileCombined = new TFile((outPath + combinedPlotsFile).c_str(),
			(Option_t*) "RECREATE");
	general_mapCombined["rate"] = new TH2F("rate", ";V_Drift ;V_Amp",
			numberOfXBins, firstXBinValue, lastXBinValue,
			(ampEnd - ampStart) / ampSteps + 1, ampStart - 0.5 * ampSteps,
			ampEnd + 0.5 * ampSteps);
	general_mapCombined["saturatedX"] = new TH2F("saturatedX", ";VDrift ;VAmp",
			numberOfXBins, firstXBinValue, lastXBinValue,
			(ampEnd - ampStart) / ampSteps + 1, ampStart - 0.5 * ampSteps,
			ampEnd + 0.5 * ampSteps);
	general_mapCombined["saturatedY"] = new TH2F("saturatedY", ";VDrift ;VAmp",
			numberOfXBins, firstXBinValue, lastXBinValue,
			(ampEnd - ampStart) / ampSteps + 1, ampStart - 0.5 * ampSteps,
			ampEnd + 0.5 * ampSteps);
	general_mapCombined["chargeX"] = new TH2F("chargeX", ";VDrift ;VAmp",
			numberOfXBins, firstXBinValue, lastXBinValue,
			(ampEnd - ampStart) / ampSteps + 1, ampStart - 0.5 * ampSteps,
			ampEnd + 0.5 * ampSteps);
	general_mapCombined["chargeY"] = new TH2F("chargeY", ";VDrift ;VAmp",
			numberOfXBins, firstXBinValue, lastXBinValue,
			(ampEnd - ampStart) / ampSteps + 1, ampStart - 0.5 * ampSteps,
			ampEnd + 0.5 * ampSteps);

// iterate of different runs in the map
	for (map<string, TFile*>::const_iterator Fitr(mapFile.begin());
			Fitr != mapFile.end(); ++Fitr) {

		int eventNumber = 0; //initialisation of counting variable for later use

		// Initializing Global Histograms
		general_mapHist1D["mmhitx"] = new TH1F("mmhitx", ";x [strips]; entries",
				xStrips, 0, xStrips);
		general_mapHist1D["mmhity"] = new TH1F("mmhity", ";y [strips]; entries",
				yStrips, 0, yStrips);
		general_mapHist1D["mmclusterx"] = new TH1F("mmclusterx",
				";x cluster size [strips]; entries", 50, 0, 50.);
		general_mapHist1D["mmclustery"] = new TH1F("mmclustery",
				";y cluster size [strips]; entries", 50, 0, 50.);
		general_mapHist1D["mmchargex"] = new TH1F("mmchargex",
				";charge X; entries", 100, 0, 1000);
		general_mapHist1D["mmchargey"] = new TH1F("mmchargey",
				";charge Y; entries", 100, 0, 1000);
		general_mapHist1D["mmtimex"] = new TH1F("mmtimex",
				";time [ns]; entries", (TRGBURST + 1) * 3, 0,
				(TRGBURST + 1) * 3 * 25.);
		general_mapHist1D["mmtimey"] = new TH1F("mmtimey",
				";time [ns]; entries", (TRGBURST + 1) * 3, 0,
				(TRGBURST + 1) * 3 * 25.);
		general_mapHist1D["mmdtime"] = new TH1F("mmdtime",
				";#Delta time [s]; entries", 500, 0, 50.);
		general_mapHist1D["mmrate"] = new TH1F("mmrate",
				";rate/10min [Hz]; entries", 200, 0, 2.);

		general_mapHist2D["mmhitmap"] = new TH2F("mmhitmap",
				";x [strips]; y [strips]", xStrips, 0, xStrips, yStrips, 0,
				yStrips);
		general_mapHist2D["mmhitmapMaxCut"] = new TH2F("mmhitmapMaxCut",
				";x [strips]; y [strips]", xStrips, 0, xStrips, yStrips, 0,
				yStrips);
		general_mapHist2D["mmhitmapGausCut"] = new TH2F("mmhitmapGausCut",
				";x [strips]; y [strips]", xStrips * 4, 0, xStrips, yStrips * 4,
				0, yStrips);

		//initialize trees with structure defined above
		TTree* fitTree = new TTree("T", "results of gauss fit");

		fitTree->Branch("gauss", &(gauss.gaussXmean),
				"gaussXmean/D:gaussXmeanError/D:gaussXsigma/D:gaussXcharge/D:gaussXchi/D:gaussXdof/D:gaussXchiRed/D:gaussYmean/D:gaussYmeanError:gaussYsigma/D:gaussYcharge/D:gaussYchi/D:gaussYdof/D:gaussYchiRed/D:number/I");
		fitTree->Branch("maxi", &maxi.maxXmean,
				"maxXmean/I:maxXcharge:maxXcluster:maxYmean:maxYcharge:maxYcluster:number");

		general_mapTree["fits"] = fitTree;

		// initialize the output file for analysis of each run
		TFile *file0 = Fitr->second;

		// Read NUTuple and execute events
		vector<string> vec_Filenames = MicroMegas->getFileName(Fitr->first);
		m_event = new MMQuickEvent(vec_Filenames, "raw", -1); //last number indicates number of events to be analysed, -1 for all events
		m_TotalEventNumber = m_event->getEventNumber();

		// loop over all events
		while (m_event->getNextEvent()) {
			if (analyseMMEvent(m_event, eventNumber, TRGBURST) == true) {
			}
			eventNumber++;
		}

		// fill dtime + rate hist
		vector<double> ratesOverMeasurementTime;
		float lengthOfMeasurement = 0.;
		sort(eventTimes.begin(), eventTimes.end());
		float tempTime = eventTimes.at(0);
		float tempDeltaTime = 0.;
		float timePeriod = 30.; // time period for rate hist (10 s)
		int periodCount = 0;
		int beginOfTimePeriod = 0;
		double lastTime = eventTimes.at(0);
		for (int e = 1; e < eventTimes.size(); e++) { // start at 1 because first is already loaded
			float deltaTime = eventTimes.at(e) - lastTime;
			if (deltaTime < 1000.) { // to get malformed events out
				general_mapHist1D["mmdtime"]->Fill(deltaTime);
				lengthOfMeasurement += deltaTime;
				tempDeltaTime += deltaTime;
				if (tempDeltaTime >= timePeriod) {
					ratesOverMeasurementTime.push_back(e - beginOfTimePeriod);
					general_mapHist1D["mmrate"]->Fill(
							(e - beginOfTimePeriod) / tempDeltaTime);
					periodCount++;
					beginOfTimePeriod = e;
					tempDeltaTime = 0.;
				}
			}
			lastTime = eventTimes.at(e);
		}
		eventTimes.clear(); // clear vector for next measurment

		//delete m_event to clear cache
		delete m_event;

//TODO: Fill histograms which show results of all measurement
//always remove "1" after the comment
		fileCombined->cd();
		general_mapCombined["rate"]->SetBinContent(
				(atoi(Fitr->first.substr(2, 3).c_str()) - driftStart)
						/ driftSteps + 1,
				(atoi(Fitr->first.substr(7, 3).c_str()) - ampStart) / ampSteps
						+ 1,/*insert here number of events*/
				1 / lengthOfMeasurement);
		general_mapCombined["saturatedX"]->SetBinContent(
				(atoi(Fitr->first.substr(2, 3).c_str()) - driftStart)
						/ driftSteps + 1,
				(atoi(Fitr->first.substr(7, 3).c_str()) - ampStart) / ampSteps
						+ 1,/*insert here number of hits with saturation in X*/
				1 / lengthOfMeasurement);
		general_mapCombined["saturatedY"]->SetBinContent(
				(atoi(Fitr->first.substr(2, 3).c_str()) - driftStart)
						/ driftSteps + 1,
				(atoi(Fitr->first.substr(7, 3).c_str()) - ampStart) / ampSteps
						+ 1,/*insert here number of hits with saturation in Y*/
				1 / lengthOfMeasurement);
		general_mapCombined["chargeX"]->SetBinContent(
				(atoi(Fitr->first.substr(2, 3).c_str()) - driftStart)
						/ driftSteps + 1,
				(atoi(Fitr->first.substr(7, 3).c_str()) - ampStart) / ampSteps
						+ 1,/*insert charge of X here*/
				1);
		general_mapCombined["chargeY"]->SetBinContent(
				(atoi(Fitr->first.substr(2, 3).c_str()) - driftStart)
						/ driftSteps + 1,
				(atoi(Fitr->first.substr(7, 3).c_str()) - ampStart) / ampSteps
						+ 1,/*insert charge of Y here*/
				1);

		/// Saving Results
		file0->cd("trees");

		/// loop over map of the plots for saving
		for (map<string, TH1F*>::iterator iter = general_mapHist1D.begin();
				iter != general_mapHist1D.end(); iter++) {
			iter->second->SetOption("error");
			iter->second->Write();
			delete iter->second;
		}
		for (map<string, TH2F*>::iterator iter = general_mapHist2D.begin();
				iter != general_mapHist2D.end(); iter++) {
			iter->second->SetOption("error");
			iter->second->Write();
			delete iter->second;
		}
		gDirectory->mkdir("2D_Events");
		gDirectory->cd("2D_Events");
		for (map<string, TH2F*>::iterator iter = general_mapHist2DEvent.begin();
				iter != general_mapHist2DEvent.end(); iter++) {
			iter->second->SetOption("error");
			iter->second->Write();
			delete iter->second;
		}
		gDirectory->cd("..");
		gDirectory->mkdir("Fits");
		gDirectory->cd("Fits");
		for (map<string, TH1F*>::iterator iter = general_mapPlotFit.begin();
				iter != general_mapPlotFit.end(); iter++) {
			iter->second->SetName(iter->first.c_str());
			iter->second->Write();
			delete iter->second;
		}
		gDirectory->cd("..");
		gDirectory->mkdir("Trees");
		gDirectory->cd("Trees");
		for (map<string, TTree*>::iterator iter = general_mapTree.begin();
				iter != general_mapTree.end(); iter++) {
			iter->second->Write();
			delete iter->second;
		}
		file0->Close();

		//end of processing of one run
	}

//save combined plots
	fileCombined->cd();
	for (map<string, TH2F*>::iterator iter = general_mapCombined.begin();
			iter != general_mapCombined.end(); iter++) {
		iter->second->Write();
		delete iter->second;
	}
	fileCombined->Close();

	return 0;
}
