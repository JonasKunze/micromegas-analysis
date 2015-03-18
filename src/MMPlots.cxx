#include "CCommonIncludes.h"
#include "MMQuickEvent.h"
#include "MapFile.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "Helper.h"

#include <thread>
#include <set>
/*
 * Limit the number of events to be processed to gain speed for debugging
 * -1 means all events will be processed
 */
int MAX_NUM_OF_EVENTS_TO_BE_PROCESSED = -1; // 192154 (run with fewest events)
#define MAX_NUM_OF_RUNS_TO_BE_PROCESSED -1

#define DRAW_CUT_EVENT_DISPLAYS true
/*
 * Cuts
 */
// Minimal charge required for the strip with maximum charge
#define MIN_CHARGE_X 60
#define MIN_CHARGE_Y 120

//#define MIN_CLUSTER_X 2
//#define MAX_CLUSTER_X 8
//#define MIN_CLUSTER_Y 2
//#define MAX_CLUSTER_Y 25

//#define MIN_CHARGE_X 0
//#define MIN_CHARGE_Y 0

#define MIN_TIMESLICE 1
#define MAX_TIMESLICE 25
#define MAX_XY_TIME_DIFFERENCE 1
#define MIN_XY_TIME_DIFFERENCE 0

#define FIT_RANGE 20
#define MAX_FIT_MEAN_DISTANCE_TO_MAX 2 // Number of strips

using namespace std;

// Global Variables
MMQuickEvent *m_event;
map<string, TTree*> general_mapTree; 		//TTress
map<string, TH1F*> general_mapHist1D; 	//1D histogram of analysis for each run
map<string, TH2F*> general_mapHist2D;	//2D histogram of analysis for each run
map<string, TH1F*> general_mapPlotFit;		//plot of fits
map<string, TH2F*> general_mapCombined;		//combined Plots
map<string, TH1F*> general_mapCombined1D;

map<string, TH2F*> global_mapCombined2D;

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

CutStatistic nocut_EventsWithSmallCharge("nocut_smallChargeEvents");
CutStatistic nocut_xtimeCutLargeYTimeEvents("nocut_xtimeCutLargeYTimeEvents");

CutStatistic timingCuts("a_timingCuts");
CutStatistic timeCoincidenceCuts("b_timeCoincidenceCuts");
CutStatistic chargeCuts("c_chargeCuts");
//CutStatistic clusterCuts("bb_clusterCuts");
CutStatistic absolutePositionXCuts("d_absolutePositionXCuts");
CutStatistic absolutePositionYCuts("e_absolutePositionYCuts");
CutStatistic proportionXCuts("f_proportionXCuts");
CutStatistic proportionYCuts("g_proportionYCuts");
CutStatistic fitProblemCuts("h_fitProblemCuts");
CutStatistic fitMeanMaxChargeDistanceCuts("i_fitMeanMaxChargeDistanceCuts");

/**
 * Returns true for every Nth eventNumber so that about totalStores times true is returned for any number of events
 */
bool storeHistogram(int eventNumber, int totalStores) {
	if (MAX_NUM_OF_EVENTS_TO_BE_PROCESSED < 0)
		return eventNumber % 600000 / totalStores == 0;

	return (MAX_NUM_OF_EVENTS_TO_BE_PROCESSED > totalStores
			&& eventNumber % (MAX_NUM_OF_EVENTS_TO_BE_PROCESSED / totalStores)
					== 0)
			|| (MAX_NUM_OF_EVENTS_TO_BE_PROCESSED > 0
					&& MAX_NUM_OF_EVENTS_TO_BE_PROCESSED <= totalStores);
}

// analysis of single event: characteristics of event and Gaussian fit
bool analyseMMEvent(MMQuickEvent *event, int eventNumber, int TRGBURST) {

// declaration of helping variables to more easily access event data
	vector<unsigned int> apvIDofStrip = *event->apv_id; // MMQuickEvent::isX(apvIDofStrip[i]) returns true if the i-th strip is X-layer
	vector<unsigned int> stripNumShowingSignal = *event->mm_strip; // stripNumShowingSignal[i] is absolute strip number (strips without charge are not stored anywhere)
	vector<vector<short> > chargeOfStripOfTime = *event->apv_q; // chargeOfStripOfTime[i][j] is the charge of strip i in time section j (matrix of whole event)
	vector<short> maxChargeOfStrip = *event->apv_qmax; // maxChargeOfStrip[i] is the maxmimal measured charge of strip i of all time sections
	vector<short> timeSliceOfMaxChargeOfStrip = *event->apv_tbqmax; // timeSliceOfMaxChargeOfStrip[i] is the time section of the corresponding maximum charge (see above)

	/*
	 * 2. Find maximum charge
	 */
	event->findMaxCharge();

	general_mapHist1D["mmchargexUncut"]->Fill(event->maxChargeX);
	general_mapHist1D["mmchargeyUncut"]->Fill(event->maxChargeY);

	general_mapCombined1D["chargexAllEventsUncut"]->Fill(event->maxChargeX);
	general_mapCombined1D["chargeyAllEventsUncut"]->Fill(event->maxChargeY);

	general_mapCombined1D["timeDistributionUncutX"]->Fill(
			event->timeSliceOfMaxChargeX);
	general_mapCombined1D["timeDistributionUncutY"]->Fill(
			event->timeSliceOfMaxChargeY);

	if (event->stripWithMaxChargeX != -1 && event->stripWithMaxChargeY != -1
			&& storeHistogram(eventNumber, 10000)) {
		event->generateTimeShape(general_mapCombined["timeShapeXUncut"],
				event->maxChargeX, event->stripWithMaxChargeX,
				event->timeSliceOfMaxChargeX);
		event->generateTimeShape(general_mapCombined["timeShapeYUncut"],
				event->maxChargeY, event->stripWithMaxChargeY,
				event->timeSliceOfMaxChargeY);
	}

	// Timing cut
	if (event->timeSliceOfMaxChargeX < MIN_TIMESLICE
			|| event->timeSliceOfMaxChargeX > MAX_TIMESLICE) {
		timingCuts.Fill(1, event);
		if (event->timeSliceOfMaxChargeX != -1
				&& event->timeSliceOfMaxChargeY > 0
				&& event->timeSliceOfMaxChargeY < 7) {
			nocut_xtimeCutLargeYTimeEvents.Fill(0, event);
		}
		return false;
	}

	general_mapCombined1D["timeDistributionYAfterTimeXCut"]->Fill(
			event->timeSliceOfMaxChargeY);

	if (event->timeSliceOfMaxChargeY < MIN_TIMESLICE
			|| event->timeSliceOfMaxChargeY > MAX_TIMESLICE) {
		timingCuts.Fill(1, event);
		return false;
	} else {
		timingCuts.Fill(0, event);
	}

	general_mapCombined1D["timeDistributionXAfterTimeCut"]->Fill(
			event->timeSliceOfMaxChargeX);
	general_mapCombined1D["timeDistributionYAfterTimeCut"]->Fill(
			event->timeSliceOfMaxChargeY);

	general_mapCombined1D["chargexAllEventsAfterTimingCut"]->Fill(
			event->maxChargeX);
	general_mapCombined1D["chargeyAllEventsAfterTimingCut"]->Fill(
			event->maxChargeY);

	if (event->timeSliceOfMaxChargeX != -1
			&& event->timeSliceOfMaxChargeY != -1) {
		general_mapCombined1D["timeCoincidence"]->Fill(
				event->timeSliceOfMaxChargeX - event->timeSliceOfMaxChargeY);
	}

	// coincidence cut
	if (event->timeSliceOfMaxChargeX
			- event->timeSliceOfMaxChargeY> MAX_XY_TIME_DIFFERENCE ||
			event->timeSliceOfMaxChargeX - event->timeSliceOfMaxChargeY
			< MIN_XY_TIME_DIFFERENCE) {

		if (event->maxChargeX < MIN_CHARGE_X || event->maxChargeY < MIN_CHARGE_Y) {
			nocut_EventsWithSmallCharge.Fill(0, event);
		}

		timeCoincidenceCuts.Fill(1, event);
		return false;
	} else {
		timeCoincidenceCuts.Fill(0, event);
	}
	general_mapCombined1D["chargexAllEventsAfterCoincidenceCut"]->Fill(
			event->maxChargeX);
	general_mapCombined1D["chargeyAllEventsAfterCoincidenceCut"]->Fill(
			event->maxChargeY);

	// Charge cut
	if (event->maxChargeX < MIN_CHARGE_X || event->maxChargeY < MIN_CHARGE_Y) {
		chargeCuts.Fill(1, event);
		return false;
	} else {
		chargeCuts.Fill(0, event);
	}

	/*
	 * Calculate cluster sizes
	 */
//	event->generateFixedTimeCrossSections();
//
//	int clusterSizeX = event->calculateClusterSize(
//			event->stripAndChargeAtMaxChargeTimeX,
//			event->positionOfMaxChargeInCrossSectionX);
//
//	int clusterSizeY = event->calculateClusterSize(
//			event->stripAndChargeAtMaxChargeTimeY,
//			event->positionOfMaxChargeInCrossSectionY);
//
//	general_mapHist1D["mmclusterxUncut"]->Fill(clusterSizeX);
//	general_mapHist1D["mmclusteryUncut"]->Fill(clusterSizeY);
//
//	general_mapCombined1D["clusterxUncut"]->Fill(clusterSizeX);
//	general_mapCombined1D["clusteryUncut"]->Fill(clusterSizeY);
//
//	// Cluster cut
//	if (clusterSizeX < MIN_CLUSTER_X || clusterSizeY < MIN_CLUSTER_Y
//			|| clusterSizeX > MAX_CLUSTER_X || clusterSizeY > MAX_CLUSTER_Y) {
//		std::stringstream suffix;
//		suffix << clusterSizeX << "-" << clusterSizeY;
//		clusterCuts.Fill(1, event, suffix.str());
//		return false;
//	} else {
//		clusterCuts.Fill(0, event);
//	}
	/*
	 * 4. Gaussian fits to charge distribution over strips at timestep with maximum charge
	 */

	event->generateFixedTimeCrossSections();
	// Proportion cuts
	bool acceptEventX = event->runProportionCut(
			general_mapCombined["mmhitneighboursX"],
			event->stripAndChargeAtMaxChargeTimeX, event->maxChargeX,
			MapFile::getProportionLimitsOfMaxHitNeighboursX(),
			absolutePositionXCuts, proportionXCuts, false,
			event->positionOfMaxChargeInCrossSectionX);

	bool acceptEventY = event->runProportionCut(
			general_mapCombined["mmhitneighboursY"],
			event->stripAndChargeAtMaxChargeTimeY, event->maxChargeY,
			MapFile::getProportionLimitsOfMaxHitNeighboursY(),
			absolutePositionYCuts, proportionYCuts, !acceptEventX,
			event->positionOfMaxChargeInCrossSectionY);

	if (!acceptEventX || !acceptEventY) {
		return false;
	}

	/*
	 * Fit hits
	 */
	TF1* gaussFitX = NULL;
	TF1* gaussFitY = NULL;
	TH1F* fitHistoX = NULL;
	TH1F* fitHistoY = NULL;

	int startFitRange = stripNumShowingSignal[event->stripWithMaxChargeX]
			- FIT_RANGE / 2;
	gaussFitX = fitGauss(event->stripAndChargeAtMaxChargeTimeX, eventNumber,
			"maxChargeCrossSectionX", fitHistoX,
			startFitRange > 0 ? startFitRange : 0,
			stripNumShowingSignal[event->stripWithMaxChargeX] + FIT_RANGE / 2);

	// fit problem cut
	if (gaussFitX == NULL) {
		fitProblemCuts.Fill(1, event);
		delete fitHistoX;
		return false;
	}
	/*
	 * Check if the fit mean is close enough to the maximum
	 */
	// fit mean distance cut
	double mean = gaussFitX->GetParameter(1);
	if (abs(
			stripNumShowingSignal[event->stripWithMaxChargeX]
					- mean) > MAX_FIT_MEAN_DISTANCE_TO_MAX) {
		delete fitHistoX;
		fitProblemCuts.Fill(0, event);
		fitMeanMaxChargeDistanceCuts.Fill(1, event);
		return false;
	}

	gaussFitY = fitGauss(event->stripAndChargeAtMaxChargeTimeY, eventNumber,
			"maxChargeCrossSectionY", fitHistoY,
			stripNumShowingSignal[event->stripWithMaxChargeY] - FIT_RANGE / 2,
			stripNumShowingSignal[event->stripWithMaxChargeY] + FIT_RANGE / 2);

	if (gaussFitY == NULL) {
		fitProblemCuts.Fill(1, event);
		delete fitHistoX;
		delete fitHistoY;
		return false;
	} else {
		fitProblemCuts.Fill(0, event);
	}

	/*
	 * Check if the fit mean is close enough to the maximum
	 */
	mean = gaussFitY->GetParameter(1);
	if (abs(
			stripNumShowingSignal[event->stripWithMaxChargeY]
					- mean) > MAX_FIT_MEAN_DISTANCE_TO_MAX) {
		delete fitHistoX;
		delete fitHistoY;
		fitMeanMaxChargeDistanceCuts.Fill(1, event);
		return false;
	} else {
		fitMeanMaxChargeDistanceCuts.Fill(0, event);
	}

	/*
	 * ############################################################
	 * #################### ALL CUTS DONE HERE ####################
	 * ############################################################
	 */
	event->generateTimeShape(general_mapCombined["timeShapeX"],
			event->maxChargeX, event->stripWithMaxChargeX,
			event->timeSliceOfMaxChargeX);
	event->generateTimeShape(general_mapCombined["timeShapeY"],
			event->maxChargeY, event->stripWithMaxChargeY,
			event->timeSliceOfMaxChargeY);

	general_mapCombined1D["timeDistributionX"]->Fill(
			event->timeSliceOfMaxChargeX);
	general_mapCombined1D["timeDistributionY"]->Fill(
			event->timeSliceOfMaxChargeY);

//storage after procession
//Fill trees	(replace 1)
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

	if (storeHistogram(eventNumber, 5)) {
		general_mapPlotFit[std::string(fitHistoX->GetName())] = fitHistoX;
		general_mapPlotFit[std::string(fitHistoY->GetName())] = fitHistoY;

		std::stringstream namePrefix;
		namePrefix << "DG" << MapFile::driftGap << "-";
		writeToPdf<TH1F>(fitHistoX, "HitWidthFits", "", namePrefix.str());
		writeToPdf<TH1F>(fitHistoY, "HitWidthFits", "", namePrefix.str());
	} else {
		delete fitHistoY;
		delete fitHistoX;
	}

	general_mapHist2D["mmhitmap"]->Fill(
			/*strip with maximum charge in X*/stripNumShowingSignal[event->stripWithMaxChargeX],
			/*strip with maximum charge in Y*/stripNumShowingSignal[event->stripWithMaxChargeY]);

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

//	general_mapHist1D["mmclusterx"]->Fill(clusterSizeX);
//	general_mapHist1D["mmclustery"]->Fill(clusterSizeY);
//
//	general_mapCombined1D["clusterx"]->Fill(clusterSizeX);
//	general_mapCombined1D["clustery"]->Fill(clusterSizeY);

	general_mapHist1D["mmtimex"]->Fill(
	/*time of maximum charge x*/event->timeSliceOfMaxChargeX * 25);
	general_mapHist1D["mmtimey"]->Fill(
	/*time of maximum charge y*/event->timeSliceOfMaxChargeY * 25);

	eventTimes.push_back(
			(double) event->time_s + (double) event->time_us / 1e6);
	return true;
}

// Main Program
void readFiles(MapFile MicroMegas, std::vector<double>& averageHitwidthsX,
		std::vector<double>& averageHitwidthsY,
		std::vector<double>& averageHitwidthsXError,
		std::vector<double>& averageHitwidthsYError,
		std::map<double/*ED*/,
				std::map<int/*VA*/,
						std::map<double/*DG*/,
								std::pair<double/*HitWIDTHs*/, double/*Error*/>>>>& hitwidthsByEdbyVaByDgX,
		std::map<double/*ED*/,
				std::map<int/*VA*/,
						std::map<double/*DG*/,
								std::pair<double/*HitWIDTHs*/, double/*Error*/>>>>& hitwidthsByEdbyVaByDgY) {


	for (auto& cutStat : CutStatistic::instances) {
		cutStat->reset();
	}

// map files to read different run of data in a row
// get data file name from MapFile.h
	map<string, TFile*> mapFile = MicroMegas.getFile();

//TRGBURST gives number of recorded timesteps (variable from data aquisition)
//timesteps = (TRGBURST+1)*3
	const int TRGBURST = 8;

//initialize file and histograms for combined output of all runs
	int numberOfXBins = (MicroMegas.driftEnd - MicroMegas.driftStart)
			/ MicroMegas.driftSteps + 1;
	double firstXBinValue = MicroMegas.driftStart - 0.5 * MicroMegas.driftSteps;
	double lastXBinValue = MicroMegas.driftEnd + 0.5 * MicroMegas.driftSteps;

	std::stringstream combinedFileName;
	combinedFileName << outPath << MicroMegas.getDriftGap()
			<< combinedPlotsFile;

	TFile* fileCombined = new TFile(combinedFileName.str().c_str(),
			(Option_t*) "RECREATE");
	general_mapCombined["rate"] = new TH2F("rate",
			";VDrift [V];VAmp [V];rate [Hz]", numberOfXBins, firstXBinValue,
			lastXBinValue,
			(MicroMegas.ampEnd - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,
			MicroMegas.ampStart - 0.5 * MicroMegas.ampSteps,
			MicroMegas.ampEnd + 0.5 * MicroMegas.ampSteps);
	general_mapCombined["chargeX"] = new TH2F("chargeX",
			";VDrift [V];VAmp [V];Charge", numberOfXBins, firstXBinValue,
			lastXBinValue,
			(MicroMegas.ampEnd - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,
			MicroMegas.ampStart - 0.5 * MicroMegas.ampSteps,
			MicroMegas.ampEnd + 0.5 * MicroMegas.ampSteps);
	general_mapCombined["chargeY"] = new TH2F("chargeY",
			";VDrift [V];VAmp [V];Charge", numberOfXBins, firstXBinValue,
			lastXBinValue,
			(MicroMegas.ampEnd - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,
			MicroMegas.ampStart - 0.5 * MicroMegas.ampSteps,
			MicroMegas.ampEnd + 0.5 * MicroMegas.ampSteps);

	general_mapCombined["chargeXfieldStrength"] = new TH2F(
			"chargeXfieldStrength",
			";drift field strength [kV/m] ;VAmp [V];charge", numberOfXBins,
			firstXBinValue / MicroMegas.driftGap,
			lastXBinValue / MicroMegas.driftGap,
			(MicroMegas.ampEnd - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,
			MicroMegas.ampStart - 0.5 * MicroMegas.ampSteps,
			MicroMegas.ampEnd + 0.5 * MicroMegas.ampSteps);
	general_mapCombined["chargeYfieldStrength"] = new TH2F(
			"chargeYfieldStrength",
			";drift field strength [kV/m] ;VAmp [V];charge", numberOfXBins,
			firstXBinValue / MicroMegas.driftGap,
			lastXBinValue / MicroMegas.driftGap,
			(MicroMegas.ampEnd - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,
			MicroMegas.ampStart - 0.5 * MicroMegas.ampSteps,
			MicroMegas.ampEnd + 0.5 * MicroMegas.ampSteps);

	general_mapCombined["chargeXuncut"] = new TH2F("chargeXuncut",
			";VDrift [V] ;VAmp [V];charge", numberOfXBins, firstXBinValue,
			lastXBinValue,
			(MicroMegas.ampEnd - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,
			MicroMegas.ampStart - 0.5 * MicroMegas.ampSteps,
			MicroMegas.ampEnd + 0.5 * MicroMegas.ampSteps);
	general_mapCombined["chargeYuncut"] = new TH2F("chargeYuncut",
			";VDrift [V] ;VAmp [V];charge", numberOfXBins, firstXBinValue,
			lastXBinValue,
			(MicroMegas.ampEnd - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,
			MicroMegas.ampStart - 0.5 * MicroMegas.ampSteps,
			MicroMegas.ampEnd + 0.5 * MicroMegas.ampSteps);

	general_mapCombined["mmhitneighboursX"] = new TH2F("mmhitneighboursX",
			";distance [strips]; relative charge [% of max];count", 13, -6.5,
			6.5, 20, 0, 100);

	general_mapCombined["mmhitneighboursY"] = new TH2F("mmhitneighboursY",
			";distance [strips]; relative charge [% of max];count", 13, -6.5,
			6.5, 20, 0, 100);

	general_mapCombined["timeShapeX"] = new TH2F("timeShapeX",
			";time section [25 ns]; charge relative to maximum strip [%]",
			2 * NUMBER_OF_TIME_SLICES + 1, -NUMBER_OF_TIME_SLICES - 0.5,
			NUMBER_OF_TIME_SLICES + 0.5, 45, -34.5, 100.5);
	general_mapCombined["timeShapeY"] = new TH2F("timeShapeY",
			";time section [25 ns]; charge relative to maximum strip [%]",
			2 * NUMBER_OF_TIME_SLICES + 1, -NUMBER_OF_TIME_SLICES - 0.5,
			NUMBER_OF_TIME_SLICES + 0.5, 45, -34.5, 100.5);
	general_mapCombined["timeShapeXUncut"] = new TH2F("timeShapeXUncut",
			";time section [25 ns]; charge relative to maximum charge [%]",
			2 * NUMBER_OF_TIME_SLICES + 1, -NUMBER_OF_TIME_SLICES - 0.5,
			NUMBER_OF_TIME_SLICES + 0.5, 45, -34.5, 100.5);
	general_mapCombined["timeShapeYUncut"] = new TH2F("timeShapeYUncut",
			";time section [25 ns]; charge relative to maximum charge [%]",
			2 * NUMBER_OF_TIME_SLICES + 1, -NUMBER_OF_TIME_SLICES - 0.5,
			NUMBER_OF_TIME_SLICES + 0.5, 45, -34.5, 100.5);

	general_mapCombined1D["chargexAllEvents"] = new TH1F("chargexAllEvents",
			";charge X; entries", 100, 0, 1000);
	general_mapCombined1D["chargeyAllEvents"] = new TH1F("chargeyAllEvents",
			";charge Y; entries", 100, 0, 1000);
	general_mapCombined1D["chargexAllEventsAfterTimingCut"] = new TH1F(
			"chargexAllEventsAfterTimingCut", ";charge X; entries", 100, 0,
			1000);
	general_mapCombined1D["chargeyAllEventsAfterTimingCut"] = new TH1F(
			"chargeyAllEventsAfterTimingCut", ";charge Y; entries", 100, 0,
			1000);
	general_mapCombined1D["chargexAllEventsAfterCoincidenceCut"] = new TH1F(
			"chargexAllEventsAfterCoincidenceCut", ";charge X; entries", 100, 0,
			1000);
	general_mapCombined1D["chargeyAllEventsAfterCoincidenceCut"] = new TH1F(
			"chargeyAllEventsAfterCoincidenceCut", ";charge Y; entries", 100, 0,
			1000);
	general_mapCombined1D["chargexAllEventsUncut"] = new TH1F(
			"chargexAllEventsUncut", ";charge X; entries", 100, 0, 1000);
	general_mapCombined1D["chargeyAllEventsUncut"] = new TH1F(
			"chargeyAllEventsUncut", ";charge Y; entries", 100, 0, 1000);

	general_mapCombined1D["hitWidthX"] = new TH1F("hitWidthX",
			";sigmaRunMittel ;entries", 50, 0, 3);
	general_mapCombined1D["hitWidthY"] = new TH1F("hitWidthY",
			";sigmaRunMittel ;entries", 50, 0, 3);

	general_mapCombined1D["timeDistributionXAfterTimeCut"] = new TH1F(
			"timeDistributionXAfterTimeCut", ";time section ;entries",
			NUMBER_OF_TIME_SLICES + 1, -1.5, 26.5);
	general_mapCombined1D["timeDistributionYAfterTimeCut"] = new TH1F(
			"timeDistributionYAfterTimeCut", ";time section ;entries",
			NUMBER_OF_TIME_SLICES + 1, -1.5, 26.5);
	general_mapCombined1D["timeDistributionYAfterTimeXCut"] = new TH1F(
			"timeDistributionYAfterTimeXCut", ";time section ;entries",
			NUMBER_OF_TIME_SLICES + 1, -1.5, 26.5);

	general_mapCombined1D["timeDistributionX"] = new TH1F("timeDistributionX",
			";time section ;entries", NUMBER_OF_TIME_SLICES + 1, -1.5, 26.5);
	general_mapCombined1D["timeDistributionY"] = new TH1F("timeDistributionY",
			";time section ;entries", NUMBER_OF_TIME_SLICES + 1, -1.5, 26.5);
	general_mapCombined1D["timeDistributionUncutX"] = new TH1F(
			"timeDistributionUncutX", ";time section ;entries",
			NUMBER_OF_TIME_SLICES + 1, -1.5, 26.5);
	general_mapCombined1D["timeDistributionUncutY"] = new TH1F(
			"timeDistributionUncutY", ";time section ;entries",
			NUMBER_OF_TIME_SLICES + 1, -1.5, 26.5);

	general_mapCombined1D["timeCoincidence"] = new TH1F("timeCoincidence",
			";time x-y [25 ns] ;entries", 11, -5.5, 5.5);

//	general_mapCombined1D["clusterx"] = new TH1F("clusterx",
//			";x cluster size [strips]; entries", 30, 0, 30.);
//	general_mapCombined1D["clustery"] = new TH1F("clustery",
//			";y cluster size [strips]; entries", 30, 0, 30.);
//	general_mapCombined1D["clusterxUncut"] = new TH1F("clusterxUncut",
//			";x cluster size [strips]; entries", 30, 0, 30.);
//	general_mapCombined1D["clusteryUncut"] = new TH1F("clusteryUncut",
//			";y cluster size [strips]; entries", 30, 0, 30.);

	int numberOfRunsToProcess = mapFile.size();
	if (MAX_NUM_OF_RUNS_TO_BE_PROCESSED < numberOfRunsToProcess
			&& MAX_NUM_OF_RUNS_TO_BE_PROCESSED > 0) {
		numberOfRunsToProcess = MAX_NUM_OF_RUNS_TO_BE_PROCESSED;
	}
	/*
	 * Data for graphs to plot the fit width vs the value of VD for every run
	 */
	std::vector<double> VDsForGraphsX;
	std::vector<double> VAsForGraphsX;
	std::vector<double> hitWidthsX;
	std::vector<double> hitWidthsXErrors;
	std::vector<double> VDsForGraphsY;
	std::vector<double> VAsForGraphsY;
	std::vector<double> hitWidthsY;
	std::vector<double> hitWidthsYErrors;

// iterate of different runs in the map
	int runNumber = 0;
	for (map<string, TFile*>::const_iterator Fitr(mapFile.begin());
			Fitr != mapFile.end(); ++Fitr) {

		if (runNumber == MAX_NUM_OF_RUNS_TO_BE_PROCESSED) {
			break;
		}

		std::cout << "Reading File " << ++runNumber << " out of "
				<< numberOfRunsToProcess << std::endl;

		int eventNumber = 0; //initialisation of counting variable for later use

		// Initializing Global Histograms
		general_mapHist1D["mmhitx"] = new TH1F("mmhitx", ";x [strips]; entries",
				xStrips, 0, xStrips);
		general_mapHist1D["mmhity"] = new TH1F("mmhity", ";y [strips]; entries",
				yStrips, 0, yStrips);
//		general_mapHist1D["mmclusterx"] = new TH1F("mmclusterx",
//				";x cluster size [strips]; entries", 50, 0, 50.);
//		general_mapHist1D["mmclustery"] = new TH1F("mmclustery",
//				";y cluster size [strips]; entries", 50, 0, 50.);
//		general_mapHist1D["mmclusterxUncut"] = new TH1F("mmclusterxUncut",
//				";x cluster size [strips]; entries", 50, 0, 50.);
//		general_mapHist1D["mmclusteryUncut"] = new TH1F("mmclusteryUncut",
//				";y cluster size [strips]; entries", 50, 0, 50.);
		general_mapHist1D["mmchargex"] = new TH1F("mmchargex",
				";charge X; entries", 100, 0, 1000);
		general_mapHist1D["mmchargey"] = new TH1F("mmchargey",
				";charge Y; entries", 100, 0, 1000);
		general_mapHist1D["mmchargexUncut"] = new TH1F("mmchargexUncut",
				";charge X; entries", 100, 0, 1000);
		general_mapHist1D["mmchargeyUncut"] = new TH1F("mmchargeyUncut",
				";charge Y; entries", 100, 0, 1000);
		general_mapHist1D["mmtimex"] = new TH1F("mmtimex",
				";time [ns]; entries", (TRGBURST + 1) * 3, 0,
				(TRGBURST + 1) * 3 * 25.);
		general_mapHist1D["mmtimey"] = new TH1F("mmtimey",
				";time [ns]; entries", (TRGBURST + 1) * 3, 0,
				(TRGBURST + 1) * 3 * 25.);
		general_mapHist1D["mmdtime"] = new TH1F("mmdtime",
				";#Delta time [s]; entries", 500, 0, 50.);
//		general_mapHist1D["mmrate"] = new TH1F("mmrate",
//				";rate/10min [Hz]; entries", 200, 0, 2.);

		general_mapHist1D["mmhitWidthX"] = new TH1F("mmhitWidthX",
				";sigma; entries", 50, 0., 3.);
		general_mapHist1D["mmhitWidthY"] = new TH1F("mmhitWidthY",
				";sigma; entries", 50, 0., 3.);

		general_mapHist2D["mmhitmap"] = new TH2F("mmhitmap",
				";x [strips]; y [strips]", xStrips, 0, xStrips, yStrips, 0,
				yStrips);

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
		int numberOfAcceptedEvents = 0;
		while (m_event->getNextEvent()
				&& eventNumber != MAX_NUM_OF_EVENTS_TO_BE_PROCESSED) {
			if (analyseMMEvent(m_event, eventNumber, TRGBURST) == true) {
				numberOfAcceptedEvents++;
			}
			eventNumber++;
		}

		/*
		 * Fit hit width histogram
		 */
		int VD = MicroMegas.getVDbyFileName(Fitr->first);
		int VA = MicroMegas.getVAbyFileName(Fitr->first);
		double VE = VD / MicroMegas.driftGap;
		TF1* hitWidthFitResultsX = fitHitWidhtHistogram(
				general_mapHist1D["mmhitWidthX"],
				general_mapCombined1D["hitWidthX"], VDsForGraphsX,
				VAsForGraphsX, hitWidthsX, hitWidthsXErrors, VD, VA);

		TF1* hitWidthFitResultsY = fitHitWidhtHistogram(
				general_mapHist1D["mmhitWidthY"],
				general_mapCombined1D["hitWidthY"], VDsForGraphsY,
				VAsForGraphsY, hitWidthsY, hitWidthsYErrors, VD, VA);

		/*
		 * Store HitWidth histograms
		 */
		std::stringstream namePrefix;
		namePrefix << "DG" << MapFile::driftGap << "-" << Fitr->first << "-";
		writeToPdf<TH1F>(general_mapHist1D["mmhitWidthX"], "HitWidthHistograms",
				"", namePrefix.str());
		writeToPdf<TH1F>(general_mapHist1D["mmhitWidthY"], "HitWidthHistograms",
				"", namePrefix.str());

		/*
		 * Store charge histograms
		 */
		namePrefix.str("");
		namePrefix << "DG" << MapFile::driftGap << "-" << Fitr->first << "-";
		writeToPdf<TH1F>(general_mapHist1D["mmchargex"], "ChargeHistograms", "",
				namePrefix.str());
		writeToPdf<TH1F>(general_mapHist1D["mmchargey"], "ChargeHistograms", "",
				namePrefix.str());
		writeToPdf<TH1F>(general_mapHist1D["mmchargexUncut"],
				"ChargeHistograms", "", namePrefix.str());
		writeToPdf<TH1F>(general_mapHist1D["mmchargeyUncut"],
				"ChargeHistograms", "", namePrefix.str());

		/*
		 * Store hitmap histogram
		 */
		namePrefix.str("");
		namePrefix << "DG" << MapFile::driftGap << "-" << Fitr->first << "-";
		writeToPdf<TH2F>(general_mapHist2D["mmhitmap"], "HitMapHistograms", "",
				namePrefix.str());

		hitwidthsByEdbyVaByDgX[(int) (VE)][VA][MapFile::driftGap] =
				std::make_pair(hitWidthFitResultsX->GetParameter(1),
						hitWidthFitResultsX->GetParError(1));

		hitwidthsByEdbyVaByDgY[(int) (VE)][VA][MapFile::driftGap] =
				std::make_pair(hitWidthFitResultsY->GetParameter(1),
						hitWidthFitResultsY->GetParError(1));

		float lengthOfMeasurement = 0.;
		if (!eventTimes.empty()) {
			// fill dtime + rate hist
			vector<double> ratesOverMeasurementTime(eventTimes.size() / 2);
			sort(eventTimes.begin(), eventTimes.end());
			float tempDeltaTime = 0.;
			float timePeriod = 30.; // time period for rate hist (10 s)
			int periodCount = 0;
			int beginOfTimePeriod = 0;
			double lastTime = eventTimes.at(0);
			for (unsigned int e = 1; e < eventTimes.size(); e++) { // start at 1 because first is already loaded
				float deltaTime = eventTimes.at(e) - lastTime;
				if (deltaTime < 1000.) { // to get malformed events out
					general_mapHist1D["mmdtime"]->Fill(deltaTime);
					lengthOfMeasurement += deltaTime;
					tempDeltaTime += deltaTime;
					if (tempDeltaTime >= timePeriod) {
						ratesOverMeasurementTime.push_back(
								e - beginOfTimePeriod);
//						general_mapHist1D["mmrate"]->Fill(
//								(e - beginOfTimePeriod) / tempDeltaTime);
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

		fileCombined->cd();
		general_mapCombined["rate"]->SetBinContent(
				(VD - MicroMegas.driftStart) / MicroMegas.driftSteps + 1,
				(VA - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,
				numberOfAcceptedEvents / lengthOfMeasurement);

		global_mapCombined2D["RateByVAED"]->Fill(VD / MicroMegas.driftGap, VA,
				numberOfAcceptedEvents / lengthOfMeasurement);

		global_mapCombined2D["hitWidthYByVAED"]->Fill(VE, VA,
				hitWidthFitResultsY->GetParameter(1));
		global_mapCombined2D["hitWidthXByVAED"]->Fill(VE, VA,
				hitWidthFitResultsX->GetParameter(1));
		global_mapCombined2D["hitWidthByVAEDCounter"]->Fill(VE, VA, 1);

		global_mapCombined2D["RateByVAEDCounter"]->Fill(
				VD / MicroMegas.driftGap, VA, 1);

		global_mapCombined2D["rateVsDriftGap"]->Fill(MicroMegas.getDriftGap(),
				numberOfAcceptedEvents / lengthOfMeasurement);

		general_mapCombined["chargeX"]->SetBinContent(
				(VD - MicroMegas.driftStart) / MicroMegas.driftSteps + 1,
				(VA - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,/*insert charge of X here*/
				general_mapHist1D["mmchargex"]->GetMean());
		general_mapCombined["chargeY"]->SetBinContent(
				(VD - MicroMegas.driftStart) / MicroMegas.driftSteps + 1,
				(VA - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,/*insert charge of Y here*/
				general_mapHist1D["mmchargey"]->GetMean());

		general_mapCombined["chargeXfieldStrength"]->SetBinContent(
				(VD - MicroMegas.driftStart) / MicroMegas.driftSteps + 1,
				(VA - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,/*insert charge of X here*/
				general_mapHist1D["mmchargex"]->GetMean());
		general_mapCombined["chargeYfieldStrength"]->SetBinContent(
				(VD - MicroMegas.driftStart) / MicroMegas.driftSteps + 1,
				(VA - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,/*insert charge of Y here*/
				general_mapHist1D["mmchargey"]->GetMean());
		general_mapCombined["chargeXuncut"]->SetBinContent(
				(VD - MicroMegas.driftStart) / MicroMegas.driftSteps + 1,
				(VA - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,/*insert charge of X here*/
				general_mapHist1D["mmchargexUncut"]->GetMean());
		general_mapCombined["chargeYuncut"]->SetBinContent(
				(VD - MicroMegas.driftStart) / MicroMegas.driftSteps + 1,
				(VA - MicroMegas.ampStart) / MicroMegas.ampSteps + 1,/*insert charge of Y here*/
				general_mapHist1D["mmchargexUncut"]->GetMean());

		/// Saving Results
		file0->cd();

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
	averageHitwidthsXError.push_back(
			general_mapCombined1D["hitWidthX"]->GetMeanError());
	averageHitwidthsY.push_back(general_mapCombined1D["hitWidthY"]->GetMean());
	averageHitwidthsYError.push_back(
			general_mapCombined1D["hitWidthY"]->GetMeanError());

	/*
	 * Print cut statistics
	 */
	std::cout << "############################################################"
			<< std::endl;
	std::cout << "########## Cut statistics for " << MicroMegas.driftGap
			<< " mm driftgap" << " ##########" << std::endl;
	std::cout << "############################################################"
			<< std::endl;
	std::cout << "Name\taccepted\tcut\t%cut" << std::endl;
	for (auto& cutStat : CutStatistic::instances) {
		TH1F histo = cutStat->counterHistogram;
		double accepted = histo.GetBinContent(1);
		double cut = histo.GetBinContent(2);
		std::cout << cutStat->getName() << "\t" << accepted << "\t" << cut
				<< "\t" << (100 * cut / (accepted + cut)) << std::endl;
	}

//save combined plots
	fileCombined->cd();
	for (map<string, TH2F*>::iterator iter = general_mapCombined.begin();
			iter != general_mapCombined.end(); iter++) {
		std::stringstream subfolder;
		subfolder << MicroMegas.driftGap;
		writeTH2FToPdf(iter->second, subfolder.str(), "colz");
		iter->second->Write();
		delete iter->second;
	}

	if (DRAW_CUT_EVENT_DISPLAYS) {
		fileCombined->cd();
		fileCombined->mkdir("Cuts");
		fileCombined->cd("Cuts");
		for (auto& cutStat : CutStatistic::instances) {
			fileCombined->cd("Cuts");
			cutStat->counterHistogram.Write();

			std::stringstream eventDisplayDirName;
			eventDisplayDirName << cutStat->getName() << "events";

			std::stringstream subdir;
			subdir << MicroMegas.driftGap << "/" << eventDisplayDirName.str()
					<< "/";

			gDirectory->mkdir(eventDisplayDirName.str().c_str());
			gDirectory->cd(eventDisplayDirName.str().c_str());
			{

				gDirectory->mkdir("Cut");
				gDirectory->cd("Cut");
				for (auto& display : cutStat->eventDisplaysCut) {
					display->Write();
					writeTH2FToPdf(display, subdir.str() + "cut", "colz");
					delete display;
				}
				cutStat->eventDisplaysCut.clear();
				gDirectory->cd("..");
			}
			{
				gDirectory->mkdir("Accepted");
				gDirectory->cd("Accepted");
				for (auto& display : cutStat->eventDisplaysAccepted) {
					display->Write();
					writeTH2FToPdf(display, subdir.str() + "accepted", "colz");
					delete display;
				}
				cutStat->eventDisplaysAccepted.clear();
				gDirectory->cd("..");
			}
			gDirectory->cd("..");
		}
	}
	fileCombined->cd("");

	/*
	 * Write any TH1F to the file and to pdf
	 */
	for (map<string, TH1F*>::iterator iter = general_mapCombined1D.begin();
			iter != general_mapCombined1D.end(); iter++) {
		iter->second->Write();

		/*
		 * write PDF
		 */
		stringstream subfolder;
		subfolder << MicroMegas.driftGap;
		writeToPdf<TH1F>(iter->second, subfolder.str(), "colz");
		delete iter->second;
	}

	/*
	 * Fit and write the graphs
	 */
	fileCombined->mkdir("Graphs");
	fileCombined->cd("Graphs");
	std::set<int> allVAs, allVDs; // TreeSet to make every entry stored only once
	allVAs.insert(VAsForGraphsX.begin(), VAsForGraphsX.end());
	allVDs.insert(VDsForGraphsX.begin(), VDsForGraphsX.end());

	for (int VA : allVAs) {
		std::stringstream name;
		name << "hitWidthVsVDX-VA" << VA;
		plotHitWidthGraph(name.str(), "VD [V]", VDsForGraphsX, hitWidthsX,
				hitWidthsXErrors, VAsForGraphsX, VA, MicroMegas.driftGap,
				20 * MicroMegas.driftGap, 2000); // Skip first bin as it's bad!
		name.str("");
		name << "hitWidthVsVDY-VA" << VA;
		plotHitWidthGraph(name.str(), "VD [V]", VDsForGraphsY, hitWidthsY,
				hitWidthsYErrors, VAsForGraphsY, VA, MicroMegas.driftGap,
				20 * MicroMegas.driftGap, 2000); // Skip first bin as it's bad!
	}

	for (int VD : allVDs) {
		std::stringstream name;
		name << "hitWidthVsVAX-VD" << VD;
		plotHitWidthGraph(name.str(), "VA [V]", VAsForGraphsX, hitWidthsX,
				hitWidthsXErrors, VDsForGraphsX, VD, MicroMegas.driftGap, 0,
				1000);
		name.str("");
		name << "hitWidthVsVAY-VD" << VD;
		plotHitWidthGraph(name.str(), "VA [V]", VAsForGraphsY, hitWidthsY,
				hitWidthsYErrors, VDsForGraphsY, VD, MicroMegas.driftGap, 0,
				1000);
	}

	fileCombined->Close();
}

/*
 * Devides every bin content in histo by the according bin content in histoCounter
 */
void calculateAveragesFromTH2F(TH2F* histo, TH2F* histoCounter) {
	for (int i = 0; i < histo->GetSize(); i++) {
		if (histoCounter->GetBinContent(i) > 0) {
			double average = histo->GetBinContent(i)
					/ histoCounter->GetBinContent(i);
			histo->SetBinContent(i, average);
		}
	}
}

void initialize() {
	int driftStart = 50;
	int driftEnd = 350;
	int driftSteps = 75;
	double driftGap = 4.5;

	int numberOfXBins = (driftEnd - driftStart) / driftSteps + 1;
	double firstXBinValue = (driftStart - 0.5 * driftSteps) / driftGap;
	double lastXBinValue = (driftEnd + 0.5 * driftSteps) / driftGap;

	global_mapCombined2D["rateVsDriftGap"] = new TH2F("rateVsDriftGap",
			";Drift gap [mm]; Rate [Hz]; Counts", 24, 4, 16, 20, 0, 500);

	global_mapCombined2D["hitWidthYByVAED"] = new TH2F("hitWidthYByVAED",
			";VD [V]; VA [V];Hit width [strips]", numberOfXBins, firstXBinValue,
			lastXBinValue, 3, 488, 563);
	global_mapCombined2D["hitWidthXByVAED"] = new TH2F("hitWidthXByVAED",
			";VD [V]; VA [V];Hit width [strips]", numberOfXBins, firstXBinValue,
			lastXBinValue, 3, 488, 563);
	global_mapCombined2D["hitWidthByVAEDCounter"] = new TH2F(
			"hitWidthByVAEDCounter", ";ED [kV/m]; VA [V]; Counts",
			numberOfXBins, firstXBinValue, lastXBinValue, 3, 488, 563);

	global_mapCombined2D["RateByVAED"] = new TH2F("RateByVAED",
			";ED [kV/m]; VA [V];Rate [Hz]", numberOfXBins, firstXBinValue,
			lastXBinValue, 3, 488, 563);
	global_mapCombined2D["RateByVAEDCounter"] = new TH2F("RateByVAEDCounter",
			";ED [kV/m]; VA [V];Counts", numberOfXBins, firstXBinValue,
			lastXBinValue, 3, 488, 563);
}
// Main Program
int main(int argc, char *argv[]) {
	if (MAX_NUM_OF_EVENTS_TO_BE_PROCESSED == -1) {
		MAX_NUM_OF_EVENTS_TO_BE_PROCESSED = 1E6; // Reduce memory consumption
	}

	std::map<std::string, std::map<int, int>> m;

	std::stringstream mkdir;
	mkdir << "mkdir -p " << outPath;
	system(mkdir.str().c_str());

	initialize();

	/*
	 * Set default style options
	 */

	std::vector<double> averageHitwidthsX;
	std::vector<double> averageHitwidthsXError;
	std::vector<double> averageHitwidthsY;
	std::vector<double> averageHitwidthsYError;

	std::map<double/*ED*/,
			std::map<int/*VA*/,
					std::map<double/*DG*/,
							std::pair<double/*HitWIDTHs*/, double/*Error*/>>>> hitwidthsByDggyVaByEdX;

	std::map<double/*ED*/,
			std::map<int/*VA*/,
					std::map<double/*DG*/,
							std::pair<double/*HitWIDTHs*/, double/*Error*/>>>> hitwidthsByDggyVaByEdY;

	std::vector<double> driftGaps = MapFile::getAvailableDriftGaps();

	/*
	 * Run over all days (drift gaps)
	 */
	for (auto& driftGap : driftGaps) {
		MapFile MicroMegas(inPath, outPath, appendName, driftGap);
		readFiles(MicroMegas, averageHitwidthsX, averageHitwidthsY,
				averageHitwidthsXError, averageHitwidthsYError,
				hitwidthsByDggyVaByEdX, hitwidthsByDggyVaByEdY);
	}

	/*
	 * Calculate average values of hitWidthYByVAVD
	 */
	calculateAveragesFromTH2F(global_mapCombined2D["hitWidthYByVAED"],
			global_mapCombined2D["hitWidthByVAEDCounter"]);
	calculateAveragesFromTH2F(global_mapCombined2D["hitWidthXByVAED"],
			global_mapCombined2D["hitWidthByVAEDCounter"]);
	calculateAveragesFromTH2F(global_mapCombined2D["RateByVAED"],
			global_mapCombined2D["RateByVAEDCounter"]);

// Set all driftgap errors to 0.1 mm
	std::vector<double> driftGapErrors;
	for (unsigned int i = 0; i < driftGaps.size(); i++) {
		driftGapErrors.push_back(0.1);
	}

	std::stringstream resultFileName;
	resultFileName << outPath << "result.root";

	TFile* fileCombined = new TFile(resultFileName.str().c_str(),
			(Option_t*) "RECREATE");

	fileCombined->cd();

	plotGraph("hitWidthVsDriftGapX", "DriftGap [mm]", driftGaps, 0.1,
			averageHitwidthsX, averageHitwidthsXError, "results", 0, 100, true);

	plotGraph("hitWidthVsDriftGapY", "DriftGap [mm]", driftGaps, 0.1,
			averageHitwidthsY, averageHitwidthsYError, "results", 0, 100, true);

	/*
	 * Plot HitWidth graphs for constant EDs
	 */
	fileCombined->cd();
	generateHitWidthVsDriftGap("hitWidthVsDriftGapX", "X",
			hitwidthsByDggyVaByEdX);
	generateHitWidthVsDriftGap("hitWidthVsDriftGapY", "Y",
			hitwidthsByDggyVaByEdY);

	fileCombined->cd();
	gStyle->SetOptStat(0);
	for (auto& pair : global_mapCombined2D) {
		pair.second->SetOption("error");
		pair.second->Write();

		writeTH2FToPdf(pair.second, "results", "colz");

		delete pair.second;
	}
	fileCombined->Close();

	/*
	 * Duck run
	 */
//	initialize();
//	MapFile MicroMegas(inPath, outPath, appendName, -1);
//	readFiles(MicroMegas, averageHitwidthsX, averageHitwidthsY,
//			averageHitwidthsXError, averageHitwidthsYError,
//			hitwidthsByDggyVaByEdX, hitwidthsByDggyVaByEdY);
}
