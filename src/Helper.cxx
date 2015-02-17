/*
 * Helper.cxx
 *
 *  Created on: Feb 8, 2015
 *      Author: kunzejo
 */

#include "Helper.h"

#include <math.h>
#include <TAttMarker.h>
#include <TAxis.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TNamed.h>
#include <iostream>

using namespace std;

void writeTH2FToPdf(TH2F* object, std::string subfolder,
		std::string drawOptions) {
	gStyle->SetOptStat(0);

	if (strcmp(object->GetName(), "hitWidthYByVAED") == 0) {
		object->SetMinimum(1.4);
	}

	if (strcmp(object->GetName(), "hitWidthXByVAED") == 0) {
		object->SetMinimum(0.8);
	}

	std::stringstream pdfName;
	pdfName << outPath << subfolder << "/";

	std::stringstream mkdir;
	mkdir << "mkdir -p " << pdfName.str();
	system(mkdir.str().c_str());

	pdfName << object->GetName() << ".pdf";

	TCanvas canvas("c", "data", 200, 10, 700, 500);
	canvas.SetRightMargin(0.15);
	object->GetZaxis()->SetTitleOffset(1.2);
	object->Draw(drawOptions.c_str());
	canvas.Print(pdfName.str().c_str(), "pdf");
	gStyle->SetOptStat(1);
}

TGraph* generateGraph(std::string name, std::string xTitle,
		std::vector<double> xValues, double xError, std::vector<double> yValues,
		std::vector<double> yErrors, double fitRangeStart, double fitRangeEnd,
		bool fitQuadratic, int fitLineColor = 1) {
	double xErrors[xValues.size()];
	for (unsigned int i = 0; i < xValues.size(); i++) {
		xErrors[i] = xError;
		std::cout << yValues[i] << "\t" << yErrors[i] << std::endl;
	}

	TGraphErrors* graph = new TGraphErrors(xValues.size(), &xValues[0],
			&yValues[0], xErrors, &yErrors[0]);
	graph->SetTitle(name.c_str());
	graph->SetName(name.c_str());
	graph->GetXaxis()->SetTitle(xTitle.c_str());
	graph->GetYaxis()->SetTitle("average hit width [strips]");
	graph->SetDrawOption("AP");
	graph->SetMarkerColor(fitLineColor);
	graph->SetMarkerStyle(21);

	std::cout << "########################################" << std::endl;
	std::cout << "Fitting " << name << std::endl;
	std::cout << "########################################" << std::endl;

	/*
	 * f(x)=a*x+b
	 * b = y1-a*x1
	 * Calculate a/b by using first and last point (rough estimation only)
	 */

	double a = 0.003;
	double b = 1.3;
//	if (yValues.size() > 2) {
//		a = (yValues[2] - yValues[1]) / (xValues[2] - xValues[1]);
//		b = yValues[1] - a * xValues[1];
//	}

	if (fitQuadratic) {
		TF1 f1("f1", "pol2", fitRangeStart, fitRangeEnd);

		f1.SetParName(0, "const.");
		f1.SetParName(1, "lin.");
		f1.SetParName(2, "quad.");
		f1.SetParameter(0, b);
		f1.SetParameter(1, a);
		f1.SetParameter(2, -0.001);

		f1.SetLineColor(fitLineColor);
		graph->Fit(&f1, "Rq");
	} else {
		TF1 f1("f1", "pol1", fitRangeStart, fitRangeEnd);

		f1.SetParName(0, "const.");
		f1.SetParName(1, "lin.");
		f1.SetParameter(0, b);
		f1.SetParameter(1, a);

		f1.SetLineColor(fitLineColor);
		graph->Fit(&f1, "Rq");
	}

	return graph;
}
void plotGraph(std::string name, std::string xTitle,
		std::vector<double> xValues, double xError, std::vector<double> yValues,
		std::vector<double> yErrors, std::string subdir, double fitRangeStart,
		double fitRangeEnd, bool fitQuadratic) {
	TGraph* graph = generateGraph(name, xTitle, xValues, xError, yValues,
			yErrors, fitRangeStart, fitRangeEnd, fitQuadratic);
	graph->Write(name.c_str());

	writeToPdf<TGraph>(graph, subdir, "AP");
}

void plotHitWidthGraph(std::string name, std::string xTitle,
		std::vector<double> xValues, std::vector<double> hitWidths,
		std::vector<double> hitWidthErrors, std::vector<double> parameters,
		double parameterValue, double driftGap, double fitRangeStart,
		double fitRangeEnd) {

	gStyle->SetStatW(0.2);

	std::vector<double> xValuesFiltered;
	std::vector<double> yValuesFiltered;
	std::vector<double> yValueErrorsFiltered;
	for (unsigned int i = 0; i < xValues.size(); i++) {
		if (parameters[i] == parameterValue) {
			xValuesFiltered.push_back(xValues[i]);
			yValuesFiltered.push_back(hitWidths[i]);
			yValueErrorsFiltered.push_back(hitWidthErrors[i]);
		}
	}
	std::stringstream subfolder;
	subfolder << driftGap;

	plotGraph(name, xTitle, xValuesFiltered, 1, yValuesFiltered,
			yValueErrorsFiltered, subfolder.str(), fitRangeStart, fitRangeEnd,
			false);
}

TF1* fitHitWidhtHistogram(TH1F* mmhitWidthHisto, TH1F* combinedWidthHisto,
		std::vector<double>& VDsForGraphs, std::vector<double>& VAsForGraphs,
		std::vector<double>& hitWidthForGraphs,
		std::vector<double>& hitWidthForGraphsError, int VD, int VA) {

	// fit histrogram maxChargeCrossSection with Gaussian distribution
	mmhitWidthHisto->Fit("gaus", "Sq");
	TF1* widthHistFitResult = mmhitWidthHisto->GetFunction("gaus");
	if (widthHistFitResult) {
		combinedWidthHisto->Fill(widthHistFitResult->GetParameter(1));
		// Plot hit width vs VD
		VDsForGraphs.push_back(VD);
		VAsForGraphs.push_back(VA);
		hitWidthForGraphs.push_back(widthHistFitResult->GetParameter(1));
		hitWidthForGraphsError.push_back(widthHistFitResult->GetParError(1));
	}

	return widthHistFitResult;
}

TF1* fitGauss(vector<std::pair<int, short> > stripAndChargeAtMaxChargeTimes,
		int eventNumber, std::string name, TH1F* &maxChargeCrossSection,
		unsigned int startFitRange, unsigned int endFitRange) {
	// Generate the title of the histogram
	stringstream histoName;
	histoName.str("");
	histoName << eventNumber << name;

	// check if any hit has been passed
	if (stripAndChargeAtMaxChargeTimes.empty()) {
		return NULL;
	}

	maxChargeCrossSection = new TH1F(histoName.str().c_str(), "; strip; charge",
			endFitRange - startFitRange + 2, startFitRange, endFitRange);

// Fill the histogram
	for (unsigned int strip = 0; strip != stripAndChargeAtMaxChargeTimes.size();
			strip++) {
		if (stripAndChargeAtMaxChargeTimes[strip].first >= startFitRange
				&& stripAndChargeAtMaxChargeTimes[strip].first <= endFitRange) {
			maxChargeCrossSection->SetBinContent(
					stripAndChargeAtMaxChargeTimes[strip].first - startFitRange
							+ 1 /* Bin 0 is underflow bin => +1 */,
					stripAndChargeAtMaxChargeTimes[strip].second);
		}
	}

// fit histrogram maxChargeDistribution with Gaussian distribution
	maxChargeCrossSection->Fit("gaus", "Sq", NULL, startFitRange, endFitRange);

// return result of Gaussian fit
	return maxChargeCrossSection->GetFunction("gaus");
}

void generateHitWidthVsDriftGap(std::string title,
		std::map<double/*ED*/,
				std::map<int/*VA*/,
						std::map<double/*DG*/,
								std::pair<double/*HitWIDTHs*/, double/*Error*/>>>> hitwidthsByDggyVaByEd) {

	for (auto& EdAndVa : hitwidthsByDggyVaByEd) {
		double Ed = EdAndVa.first;

		int lineColor = 1;
		std::vector<TGraph*> graphs;
		for (auto& VaAndDg : EdAndVa.second) {
			std::vector<double> driftGaps;
			std::vector<double> HitWidths;
			std::vector<double> HitWidthErrors;

			int Va = VaAndDg.first;

			for (auto& DgAndHitWidth : VaAndDg.second) {
				driftGaps.push_back(DgAndHitWidth.first);
				HitWidths.push_back(DgAndHitWidth.second.first);
				HitWidthErrors.push_back(DgAndHitWidth.second.second);
			}

			std::stringstream title;
			title << "VA" << Va;

			std::stringstream graphSubDir;
			graphSubDir << "results/ED" << Ed;

			graphs.push_back(
					generateGraph(title.str(), "DriftGap [mm]", driftGaps, 0.1,
							HitWidths, HitWidthErrors, 0, 100, true, lineColor++));
			plotGraph(title.str(), "DriftGap [mm]", driftGaps, 0.1, HitWidths,
					HitWidthErrors, graphSubDir.str(), 0, 100, true);
		}

		/*
		 * Combine every every graph for each DG to one multigraph
		 */
		std::stringstream name;
		name << title << "-ED" << Ed;
		TMultiGraph* multigraph = new TMultiGraph(name.str().c_str(),
				";drift Gap [mm]; average hit width [strips]");
		lineColor = 1;
		for (auto& graph : graphs) {
			graph->SetMarkerColor(lineColor);
			graph->SetLineColor(lineColor);
			graph->SetMarkerStyle(19 + lineColor++);
			graph->SetFillStyle(0);
			graph->SetFillColor(0);
			multigraph->Add(graph, "");
		}

		multigraph->Write(multigraph->GetName());
		writeToPdf<TMultiGraph>(multigraph, "results", "ap", "", 0, true);
		//delete multigraph; // Why the **** can't I delete it?
	}
}
