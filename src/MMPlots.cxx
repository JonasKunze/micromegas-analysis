#include "CCommonIncludes.h"
#include "MMQuickEvent.h"
#include "MapFile.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
/*
 * Limit the number of events to be processed to gain speed for debugging
 * -1 means all events will be processed
 */
#define MAX_NUM_OF_EVENTS_TO_BE_PROCESSED 500
#define MAX_NUM_OF_RUNS_TO_BE_PROCESSED 5

/*
 * Cuts
 */
// Minimal charge required for the strip with maximum charge
#define MIN_CHARGE_X 50
#define MIN_CHARGE_Y 120
//#define MIN_CHARGE_X 0map<string, TTree*> general_mapTree; 		//TTress

//#define MIN_CHARGE_Y 0

#define MIN_TIMESLICE 3
#define MAX_TIMESLICE 23
#define MAX_XY_TIME_DIFFERENCE 1

#define FIT_RANGE 20
#define MAX_FIT_MEAN_DISTANCE_TO_MAX 2 // Number of strips

#define RUN_FITS true

using namespace std;

// Global Variables
MMQuickEvent *m_event;
map<string, TTree*> general_mapTree; 		//TTress
map<string, TH1F*> general_mapHist1D; 	//1D histogram of analysis for each run
map<string, TH2F*> general_mapHist2D;	//2D histogram of analysis for each run
map<string, TH2F*> general_mapHist2DEvent; 	//plot of 2D EventDisplay
map<string, TH1F*> general_mapPlotFit;		//plot of fits
map<string, TH2F*> general_mapCombined;		//combined Plots
map<string, TH1F*> general_mapCombined1D;
map<string, TTree*> general_mapCombinedTree;

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
cutStatistics_t cutStatistics;

//set output path and name of output files
const string inPath = "/localscratch/praktikum/data/";		//Path of the Input
const string outPath = "/localscratch/praktikum/output/"; // Path of the Output
//const string outPath = "/tmp/"; // Path of the Output
const string appendName = "";					// Name of single measurements
const string combinedPlotsFile = "combined.root";// Name of the file for the combined results of all runs (hier muss jeder Tag einzeln analysiert werden! Da Zeile 79-84(driftStart...ampSteps) für jeden Tag anders war. Es können unter anderem angeschaut werden Raten in abhängigkeit der Spannung

//Voltage range, needed for initialization of combined histograms (hier die Schritte von VD und VA angeben (amp stimmt schon)
const int driftStart = 100;
const int driftEnd = 400;
const int driftSteps = 50;
const int ampStart = 500;
const int ampEnd = 550;
const int ampSteps = 25;

/**
 * Returns true for every Nth eventNumber so that about 100-200 times true is returned for any number of events
 */
bool storeHistogram(int eventNumber) {
	return ((MAX_NUM_OF_EVENTS_TO_BE_PROCESSED < 0 && eventNumber % 5000 == 0) /* every 5000th */
			|| (MAX_NUM_OF_EVENTS_TO_BE_PROCESSED > 1000
					&& eventNumber % (MAX_NUM_OF_EVENTS_TO_BE_PROCESSED / 1000)
							== 0)
			|| (MAX_NUM_OF_EVENTS_TO_BE_PROCESSED > 0
					&& MAX_NUM_OF_EVENTS_TO_BE_PROCESSED <= 1000));
}

/**
 * Generates a new 2D histogram (heatmap) showing all measured charges in all strips of x or y direction of all times slices.
 * The histogram will be stored at general_mapHist2DEvent[eventNumber+"nameOfHistogram"]
 */
void generateEventDisplay(MMQuickEvent* event) {
	vector<vector<short> > chargeOfTimeOfStrip = *event->apv_q;
	unsigned int numberOfTimeSlices = chargeOfTimeOfStrip[0].size();

	/*
	 * Generate a new 2D histogram for the event display (x=strip, y=timeslice, z=charge)
	 */

	// Generate the title of the histogram
	stringstream histoName;
	histoName.str("");
	histoName << event->getCurrentEventNumber() << "-Eventdisplay";

	string histoNameX = histoName.str() + "_X";
	string histoNameY = histoName.str() + "_Y";

	/*
	 * Initialize a new root TH2F histogram with the right title, labels and binning:
	 *
	 */
	TH2F* eventDisplayX = new TH2F(histoNameX.c_str(),
			";Strip Number; Time [25 ns]", xStrips, 0, xStrips - 1,
			numberOfTimeSlices, 0, numberOfTimeSlices - 1);

	TH2F* eventDisplayY = new TH2F(histoNameY.c_str(),
			";Strip Number; Time [25 ns]", yStrips, 0, yStrips - 1,
			numberOfTimeSlices, 0, numberOfTimeSlices - 1);

	// Store the new histogram in the global map
	general_mapHist2DEvent[histoNameX] = eventDisplayX;
	general_mapHist2DEvent[histoNameY] = eventDisplayY;

	/*
	 * Fill eventDisplay with all measured charges at all strips for all times;
	 * for all timeSlices, we iterate through all Strips (X and Y). Depending
	 * on the apvID, the corresponding histogram is filled (X resp. Y).
	 */
	for (unsigned int timeSlice = 0; timeSlice != numberOfTimeSlices;
			timeSlice++) {
		for (unsigned int stripNum = 0; stripNum != chargeOfTimeOfStrip.size();
				stripNum++) {
			short charge = chargeOfTimeOfStrip[stripNum][timeSlice];
			unsigned int apvID = (*event->apv_id)[stripNum];
			if (MMQuickEvent::isX(apvID)) {

				// store charge of current strip in bin corresponding to absolute strip number
				eventDisplayX->SetBinContent(
						event->mm_strip->at(stripNum) + 1/*x*/,
						timeSlice + 1/*y*/, charge/*z*/);
			} else {
				eventDisplayY->SetBinContent(
						event->mm_strip->at(stripNum) + 1/*x*/,
						timeSlice + 1/*y*/, charge/*z*/);
			}
		}
	}
}

TF1* fitGauss(
		vector<std::pair<unsigned int, short> > stripAndChargeAtMaxChargeTimes,
		int eventNumber, std::string name, TH1F* &maxChargeDistribution,
		unsigned int startFitRange, unsigned int endFitRange) {
	// Generate the title of the histogram
	stringstream histoName;
	histoName.str("");
	histoName << eventNumber << name;

	// check if any hit has been passed
	if (stripAndChargeAtMaxChargeTimes.empty()) {
		return NULL;
	}

	maxChargeDistribution = new TH1F(histoName.str().c_str(), "; strip; charge",
			endFitRange - startFitRange + 2, startFitRange, endFitRange);

	// No idea why...but this needs to be done...Damn root
//	maxChargeDistribution->SetDirectory(0);
//	TH1::AddDirectory(kFALSE);

	// Fill the histogram
	for (unsigned int strip = 0; strip != stripAndChargeAtMaxChargeTimes.size();
			strip++) {
		if (stripAndChargeAtMaxChargeTimes[strip].first >= startFitRange
				&& stripAndChargeAtMaxChargeTimes[strip].first <= endFitRange) {
			maxChargeDistribution->SetBinContent(
					stripAndChargeAtMaxChargeTimes[strip].first - startFitRange
							+ 1 /* Bin 0 is underflow bin => +1 */,
					stripAndChargeAtMaxChargeTimes[strip].second);
		}
	}

// fit histrogram maxChargeDistribution with Gaussian distribution
	maxChargeDistribution->Fit("gaus", "Sq", NULL, startFitRange, endFitRange);

// return result of Gaussian fit
	return maxChargeDistribution->GetFunction("gaus");
}

// analysis of single event: characteristics of event and Gaussian fit
bool analyseMMEvent(MMQuickEvent *event, int eventNumber, int TRGBURST) {

// declaration of helping variables to more easily access event data
	vector<unsigned int> apvIDofStrip = *event->apv_id; // MMQuickEvent::isX(apvIDofStrip[i]) returns true if the i-th strip is X-layer
	vector<unsigned int> stripNumShowingSignal = *event->mm_strip; // stripNumShowingSignal[i] is absolute strip number (strips without charge are not stored anywhere)
	vector<vector<short> > chargeOfStripOfTime = *event->apv_q; // chargeOfStripOfTime[i][j] is the charge of strip i in time slice j (matrix of whole event)
	vector<short> maxChargeOfStrip = *event->apv_qmax; // maxChargeOfStrip[i] is the maxmimal measured charge of strip i of all time slices
	vector<short> timeSliceOfMaxChargeOfStrip = *event->apv_tbqmax; // timeSliceOfMaxChargeOfStrip[i] is the time slice of the corresponding maximum charge (see above)

	/*
	 * 1. Create event display
	 *
	 * Reduce the number of event display to a reasonable number
	 */
	if (storeHistogram(event->getCurrentEventNumber())) {
		generateEventDisplay(event);
	}

	/*
	 * 2. Find maximum charge
	 */
	event->findMaxCharge();

	/*
	 * zu cuts: 1. eventdisplays anschauen, erste überlegungen zu cuts
	 * 2. ladungsverteilung aller events anschauen, daran cuts festmachen (Anzahl hits gegen ladungshöhe)
	 * 3. eventuell nur ein-hit-events aussuchen, je nachdem, wie gut 2. hit von erstem zu unterscheiden ist (pulshöhendifferenz)
	 */

	/*
	 * Timing cut
	 */
	if (event->timeSliceOfMaxChargeX < MIN_TIMESLICE
			|| event->timeSliceOfMaxChargeX > MAX_TIMESLICE
			|| event->timeSliceOfMaxChargeY < MIN_TIMESLICE
			|| event->timeSliceOfMaxChargeY > MAX_TIMESLICE) {
		cutStatistics.timingCuts++;
		return false;
	}

	// coincidence check
	if (abs(
			event->timeSliceOfMaxChargeX
					- event->timeSliceOfMaxChargeY) > MAX_XY_TIME_DIFFERENCE) {
		cutStatistics.timeCoincidenceCuts++;
		return false;
	}

	/*
	 * 4. Gaussian fits to charge distribution over strips at timestep with maximum charge
	 */

// first data cut: remove events with small charge
	if (event->maxChargeX < MIN_CHARGE_X && event->maxChargeY < MIN_CHARGE_Y) {
		cutStatistics.chargeCuts++;
		return false;
	}

	event->generateFixedTimeCrossSections();

	bool acceptEventX = event->runProportionCut(
			general_mapCombined["mmhitneighboursX"],
			event->stripAndChargeAtMaxChargeTimeX, event->maxChargeX,
			MapFile::getProportionLimitsOfMaxHitNeighboursX(), cutStatistics);
	bool acceptEventY = event->runProportionCut(
			general_mapCombined["mmhitneighboursY"],
			event->stripAndChargeAtMaxChargeTimeY, event->maxChargeY,
			MapFile::getProportionLimitsOfMaxHitNeighboursY(), cutStatistics);

	if (!acceptEventX || !acceptEventY) {
		return false;
	}

	TF1* gaussFitX = NULL;
	TF1* gaussFitY = NULL;
	TH1F* fitHistoX = NULL;
	TH1F* fitHistoY = NULL;
	if (RUN_FITS) {
		int startFitRange = stripNumShowingSignal[event->stripWithMaxChargeX]
				- FIT_RANGE / 2;
		gaussFitX = fitGauss(event->stripAndChargeAtMaxChargeTimeX, eventNumber,
				"maxChargeDistributionX", fitHistoX,
				startFitRange > 0 ? startFitRange : 0,
				stripNumShowingSignal[event->stripWithMaxChargeX]
						+ FIT_RANGE / 2);

		if (gaussFitX == NULL) {
			cutStatistics.fitProblems++;
			delete fitHistoX;
			return false;
		}
		/*
		 * Check if the fit mean is close enough to the maximum
		 */
		double mean = gaussFitX->GetParameter(1);
		if (abs(
				stripNumShowingSignal[event->stripWithMaxChargeX]
						- mean) < MAX_FIT_MEAN_DISTANCE_TO_MAX) {
			if (storeHistogram(eventNumber)) {
				general_mapPlotFit[std::string(fitHistoX->GetName())] =
						fitHistoX;
			} else {
				delete fitHistoX;
				fitHistoX = NULL;
			}
		} else {
			delete fitHistoX;
			fitHistoX = NULL;
			cutStatistics.fitMeanMaxChargeDistanceCuts++;
		}

		gaussFitY = fitGauss(event->stripAndChargeAtMaxChargeTimeY, eventNumber,
				"maxChargeDistributionY", fitHistoY,
				stripNumShowingSignal[event->stripWithMaxChargeY]
						- FIT_RANGE / 2,
				stripNumShowingSignal[event->stripWithMaxChargeY]
						+ FIT_RANGE / 2);

		if (gaussFitY == NULL) {
			cutStatistics.fitProblems++;
			delete fitHistoY;
			return false;
		}

		/*
		 * Check if the fit mean is close enough to the maximum
		 */
		mean = gaussFitY->GetParameter(1);
		if (abs(
				stripNumShowingSignal[event->stripWithMaxChargeY]
						- mean) < MAX_FIT_MEAN_DISTANCE_TO_MAX) {
			if (storeHistogram(eventNumber)) {
				general_mapPlotFit[std::string(fitHistoY->GetName())] =
						fitHistoY;
			} else {
				delete fitHistoY;
				fitHistoY = NULL;
			}
		} else {
			delete fitHistoY;
			fitHistoY = NULL;
			cutStatistics.fitMeanMaxChargeDistanceCuts++;
		}
	}

//storage after procession
//Fill trees	(replace 1)
	if (/*condition to store the fit*/fitHistoY != NULL && fitHistoX != NULL) {
		gauss.gaussXmean = gaussFitX->GetParameter(1);
		gauss.gaussXmeanError = gaussFitX->GetParError(1);
		gauss.gaussXsigma = gaussFitX->GetParameter(2);
		gauss.gaussXcharge = gaussFitX->GetParameter(0);
		gauss.gaussXchi = gaussFitX->GetChisquare();
		gauss.gaussXdof = gaussFitX->GetNDF();
		gauss.gaussXchiRed = gaussFitX->GetChisquare() / gaussFitX->GetNDF();
		gauss.gaussYmean = gaussFitY->GetParameter(1);
		gauss.gaussYmeanError = gaussFitY->GetParError(1);
		gauss.gaussYsigma = gaussFitY->GetParameter(2);
		gauss.gaussYcharge = gaussFitY->GetParameter(0);
		gauss.gaussYchi = gaussFitY->GetChisquare();
		gauss.gaussYdof = gaussFitY->GetNDF();
		gauss.gaussYchiRed = gaussFitY->GetChisquare() / gaussFitY->GetNDF();
		gauss.number = eventNumber;

		general_mapHist1D["mmhitWidthX"]->Fill(gauss.gaussXsigma);
		general_mapHist1D["mmhitWidthY"]->Fill(gauss.gaussYsigma);

		/*
		 * ???
		 * Was ist hier zu tun?
		 */
		maxi.maxXmean = event->maxChargeX;
		maxi.maxYmean = event->maxChargeY;
		maxi.maxXcharge = 1;
		maxi.maxYcharge = 1;
		maxi.maxXcluster = 1;
		maxi.maxYcluster = 1;
		maxi.number = eventNumber;

		general_mapTree["fits"]->Fill();

		/*
		 * ???
		 * Was ist hier gefragt? Was ist der unterschied zu mmhitmap
		 */
		general_mapHist2D["mmhitmapMaxCut"]->Fill(/*x value*/1,/*y value*/
		1);
		general_mapHist2D["mmhitmapGausCut"]->Fill(/*x value*/1,/*y value*/
		1);

	}

	if (/*condition to fill general histograms*/!RUN_FITS
			|| (gaussFitY != NULL && gaussFitX != NULL)) {
		// coincidence check between x and y signal (within 25ns)
		general_mapHist2D["mmhitmap"]->Fill(
				/*strip with maximum charge in X*/stripNumShowingSignal[event->stripWithMaxChargeX],/*strip with maximum charge in Y*/
				stripNumShowingSignal[event->stripWithMaxChargeY]);

		general_mapCombined1D["chargexAllEvents"]->Fill(event->maxChargeX);
		general_mapCombined1D["chargeyAllEvents"]->Fill(event->maxChargeY);

		general_mapHist1D["mmchargex"]->Fill(
		/*maximum charge x*/event->maxChargeX);
		general_mapHist1D["mmchargey"]->Fill(
		/*maximum charge y*/event->maxChargeY);
		general_mapHist1D["mmhitx"]->Fill(
				/*strip x with maximum charge*/stripNumShowingSignal[event->stripWithMaxChargeX]);
		general_mapHist1D["mmhity"]->Fill(
				/*strip y with maximum charge*/stripNumShowingSignal[event->stripWithMaxChargeY]);

		/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		 * !!!!!!!!!!!!!!!!!!!!!! To be implemented!!!!!!!!!!!!!!!!!!!!!!!!!
		 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		 */
		general_mapHist1D["mmclusterx"]->Fill(
		/*number of x strips hit by one event*/event->numberOfXHits);
		general_mapHist1D["mmclustery"]->Fill(
		/*number of y strips hit by one event*/event->numberOfYHits);

		general_mapHist1D["mmtimex"]->Fill(
		/*time of maximum charge x*/event->timeSliceOfMaxChargeX * 25);
		general_mapHist1D["mmtimey"]->Fill(
		/*time of maximum charge y*/event->timeSliceOfMaxChargeY * 25);

		eventTimes.push_back(
				(double) event->time_s + (double) event->time_us / 1e6);
	}
	return true;
}

// Main Program
void readFiles(MapFile MicroMegas, std::vector<double>& averageHitwidthsX,
		std::vector<double>& averageHitwidthsY) {

// map files to read different run of data in a row
// get data file name from MapFile.h
	map<string, TFile*> mapFile = MicroMegas.getFile();

//TRGBURST gives number of recorded timesteps (variable from data aquisition)
//timesteps = (TRGBURST+1)*3
	const int TRGBURST = 8;

//initialize file and histograms for combined output of all runs
	int numberOfXBins = (driftEnd - driftStart) / driftSteps + 1;
	double firstXBinValue = driftStart - 0.5 * driftSteps;
	double lastXBinValue = driftEnd + 0.5 * driftSteps;

	std::stringstream combinedFileName;
	combinedFileName << outPath << MicroMegas.getDriftGap()
			<< combinedPlotsFile;

	std::cout << combinedFileName.str() << std::endl;
	TFile* fileCombined = new TFile(combinedFileName.str().c_str(),
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

	general_mapCombined["mmhitneighboursX"] = new TH2F("mmhitneighboursX",
			";distance [strips]; relative charge [% of max]", 21, -10, 10, 21,
			0, 100);

	general_mapCombined["mmhitneighboursY"] = new TH2F("mmhitneighboursY",
			";distance [strips]; relative charge [% of max]", 21, -10, 10, 21,
			0, 100);

	general_mapCombined1D["chargexAllEvents"] = new TH1F("chargexAllEvents",
			";charge X; entries", 100, 0, 1000);
	general_mapCombined1D["chargeyAllEvents"] = new TH1F("chargeyAllEvents",
			";charge Y; entries", 100, 0, 1000);

	general_mapCombined1D["hitWidthX"] = new TH1F("hitWidthX",
			";sigmaRunMittel ;entries", 50, 0, 3);
	general_mapCombined1D["hitWidthY"] = new TH1F("hitWidthY",
			";sigmaRunMittel ;entries", 50, 0, 3);

	//initialize trees with structure defined above
	TTree* cutTree = new TTree("Cuts", "Statistics of cuts");
	cutTree->Branch("cutStatistics", &cutStatistics,
			"timingCuts/I:chargeCuts/I:timeCoincidenceCuts/I:absolutePositionCuts/I:proportionCuts/I:fitMeanMaxChargeDistanceCuts/I:fitProblems/I");
	general_mapCombinedTree["cuts"] = cutTree;

// iterate of different runs in the map
	int runNumber = 0;
	for (map<string, TFile*>::const_iterator Fitr(mapFile.begin());
			Fitr != mapFile.end(); ++Fitr) {

		if (runNumber == MAX_NUM_OF_RUNS_TO_BE_PROCESSED - 1) {
			break;
		}

		std::cout << "Reading File " << ++runNumber << " out of "
				<< mapFile.size() << std::endl;

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

		general_mapHist1D["mmhitWidthX"] = new TH1F("mmhitWidthX",
				";sigma; entries", 30, 0., 5.);
		general_mapHist1D["mmhitWidthY"] = new TH1F("mmhitWidthY",
				";sigma; entries", 30, 0., 5.);

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
		vector<string> vec_Filenames = MicroMegas.getFileName(Fitr->first);
		m_event = new MMQuickEvent(vec_Filenames, "raw", -1); //last number indicates number of events to be analysed, -1 for all events
		m_TotalEventNumber = m_event->getEventNumber();

		// loop over all events
		while (m_event->getNextEvent()
				&& eventNumber != MAX_NUM_OF_EVENTS_TO_BE_PROCESSED) {
			cutStatistics = {0};
			if (analyseMMEvent(m_event, eventNumber, TRGBURST) == true) {
			}
			eventNumber++;
			general_mapCombinedTree["cuts"]->Fill();
		}

		/*
		 * Fit hit width histogram
		 */
		TH1F* widthHistoX = general_mapHist1D["mmhitWidthX"];
		// fit histrogram maxChargeDistribution with Gaussian distribution
		widthHistoX->Fit("gaus", "Sq");
		TF1* widthHistXFitResult = widthHistoX->GetFunction("gaus");
		if (widthHistXFitResult)
			general_mapCombined1D["hitWidthX"]->Fill(
					widthHistXFitResult->GetParameter(1));

		TH1F* widthHistoY = general_mapHist1D["mmhitWidthY"];
		// fit histrogram maYChargeDistribution with Gaussian distribution
		widthHistoY->Fit("gaus", "Sq");

		TF1* widthHistYFitResult = widthHistoY->GetFunction("gaus");
		if (widthHistXFitResult)
			general_mapCombined1D["hitWidthY"]->Fill(
					widthHistYFitResult->GetParameter(1));

		float lengthOfMeasurement = 0.;
		if (!eventTimes.empty()) {
			// fill dtime + rate hist
			vector<double> ratesOverMeasurementTime;
			sort(eventTimes.begin(), eventTimes.end());
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
						ratesOverMeasurementTime.push_back(
								e - beginOfTimePeriod);
						general_mapHist1D["mmrate"]->Fill(
								(e - beginOfTimePeriod) / tempDeltaTime);
						periodCount++;
						beginOfTimePeriod = e;
						tempDeltaTime = 0.;
					}
				}
				lastTime = eventTimes.at(e);
			}
			eventTimes.clear(); // clear vector for next measurement
		}

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
				eventNumber / lengthOfMeasurement);
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
		file0->mkdir("trees");
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
			general_mapPlotFit.erase(iter->first);
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

	/*
	 * Store the average hit widths of all runs with the current driftGap
	 */
	averageHitwidthsX.push_back(general_mapCombined1D["hitWidthX"]->GetMean());
	averageHitwidthsY.push_back(general_mapCombined1D["hitWidthY"]->GetMean());

	gDirectory->cd("..");
	gDirectory->mkdir("Combined");
	gDirectory->cd("Combined");
//save combined plots
	fileCombined->cd();
	for (map<string, TH2F*>::iterator iter = general_mapCombined.begin();
			iter != general_mapCombined.end(); iter++) {
		iter->second->Write();
		delete iter->second;
	}
	for (map<string, TH1F*>::iterator iter = general_mapCombined1D.begin();
			iter != general_mapCombined1D.end(); iter++) {
		iter->second->Write();
		delete iter->second;
	}

	gDirectory->cd("..");
	gDirectory->mkdir("Trees");
	gDirectory->cd("Trees");
	for (map<string, TTree*>::iterator iter = general_mapCombinedTree.begin();
			iter != general_mapCombinedTree.end(); iter++) {
		iter->second->Write();
		delete iter->second;
	}

	fileCombined->Close();
}

// Main Program
int main(int argc, char *argv[]) {
	std::vector<double> averageHitwidthsX;
	std::vector<double> averageHitwidthsY;
	std::vector<double> driftGaps = MapFile::getAvailableDriftGaps();

	for (std::vector<double>::iterator iter = driftGaps.begin();
			iter != driftGaps.end(); iter++) {
		MapFile MicroMegas(inPath, outPath, appendName, *iter);
		readFiles(MicroMegas, averageHitwidthsX, averageHitwidthsY);
	}

	std::stringstream resultFileName;
	resultFileName << outPath << "result.root";

	TFile* fileCombined = new TFile(resultFileName.str().c_str(),
			(Option_t*) "RECREATE");

	fileCombined->cd();
	TGraph* hitWidthVsDriftGapX = new TGraph(driftGaps.size(), &driftGaps[0],
			&averageHitwidthsX[0]);
	hitWidthVsDriftGapX->Write("hitWidthVsDriftGapX");
	TGraph* hitWidthVsDriftGapY = new TGraph(driftGaps.size(), &driftGaps[0],
			&averageHitwidthsY[0]);
	hitWidthVsDriftGapY->Write("hitWidthVsDriftGapY");
	fileCombined->Close();
}
