/*
 * Helper.cxx
 *
 *  Created on: Feb 8, 2015
 *      Author: kunzejo
 */
#ifndef HELPER_H
#define HELPER_H

#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <cstdlib>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <map>

class TF1;
class TH1;
class TMultiGraph;

class TH2F;

//set output path and name of output files
const std::string inPath = "/localscratch/praktikum/data/";	//Path of the Input
const std::string outPath = "/localscratch/praktikum/output/"; // Path of the Output
//const string outPath = "/tmp/output/"; // Path of the Output
const std::string appendName = "";				// Name of single measurements
const std::string combinedPlotsFile = "combined.root";// Name of the file for the combined results of all runs (hier muss jeder Tag einzeln analysiert werden! Da Zeile 79-84(driftStart...ampSteps) für jeden Tag anders war.

template<typename T>
void writeToPdf(T* object, std::string subfolder, std::string drawOptions,
		std::string namePrefix = "", int optFit = 1111,
		bool buildLegend = false) {
	gStyle->SetOptFit(optFit);

	std::string histoName(object->GetName());
	gStyle->SetStatW(0.2);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.47);

	if (dynamic_cast<TH1*>(object)) {
		gStyle->SetStatY(0.9);
		gStyle->SetStatX(0.9);
	}
	if (histoName == "timeDistribution") {
		gStyle->SetStatY(0.5);
		gStyle->SetStatX(0.6);
	}
	std::stringstream pdfName;
	pdfName << outPath << subfolder << "/";

	std::stringstream mkdir;
	mkdir << "mkdir -p " << pdfName.str();
	system(mkdir.str().c_str());

	pdfName << namePrefix << object->GetName() << ".pdf";

	TCanvas canvas("c", "data", 200, 10, 700, 500);
	canvas.SetLeftMargin(0.11);
	if (!dynamic_cast<TMultiGraph*>(object)) {
		object->GetYaxis()->SetTitleOffset(1.5);
	}
	object->Draw(drawOptions.c_str());
	if (buildLegend) {
		TLegend *leg = canvas.BuildLegend();
		leg->SetFillStyle(0);
		leg->SetX1(0.15);
		leg->SetY1(0.7);
		leg->SetX2(0.48);
		leg->SetY2(0.9);
	}
	canvas.Print(pdfName.str().c_str(), "pdf");
}
void writeTH2FToPdf(TH2F* object, std::string subfolder,
		std::string drawOptions);

/*
 * xValues, hitWidhts, hitWidthErrors and parameters must all have the same number of values
 *
 * (xValues[i]+-1, hitWidths[i]+-hitWidthErros[i]) will be plotted for all i with parameters[i]==parameterValue
 */
void plotHitWidthGraph(std::string name, std::string xTitle,
		std::vector<double> xValues, std::vector<double> hitWidths,
		std::vector<double> hitWidthErrors, std::vector<double> parameters,
		double parameterValue, double driftGap, double fitRangeStart,
		double fitRangeEnd);

TGraph generateGraph(std::string name, std::string xTitle,
		std::vector<double> xValues, double xError, std::vector<double> yValues,
		std::vector<double> yErrors, double fitRangeStart, double fitRangeEnd, int fitLineColor);

void plotGraph(std::string name, std::string xTitle,
		std::vector<double> xValues, double xError, std::vector<double> yValues,
		std::vector<double> yErrors, std::string subdir, double fitRangeStart,
		double fitRangeEnd);

TF1* fitHitWidhtHistogram(TH1F* mmhitWidthHisto, TH1F* combinedWidthHisto,
		std::vector<double>& VDsForGraphs, std::vector<double>& VAsForGraphs,
		std::vector<double>& hitWidthForGraphs,
		std::vector<double>& hitWidthForGraphsError, int VD, int VA);

TF1* fitGauss(
		std::vector<std::pair<unsigned int, short> > stripAndChargeAtMaxChargeTimes,
		int eventNumber, std::string name, TH1F* &maxChargeCrossSection,
		unsigned int startFitRange, unsigned int endFitRange);

void generateHitWidthVsDriftGap(std::string title,
		std::map<double/*ED*/,
				std::map<int/*VA*/,
						std::map<double/*DG*/,
								std::pair<double/*HitWIDTHs*/, double/*Error*/>>>> hitwidthsByDggyVaByEd);

#endif
