/*
 * Helper.cxx
 *
 *  Created on: Feb 8, 2015
 *      Author: kunzejo
 */

#include "Helper.h"

#include <TAttMarker.h>
#include <TAxis.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TNamed.h>
#include <iostream>
#include <vector>

void writeTH2FToPdf(TH2F* object, std::string subfolder,
		std::string drawOptions) {
	gStyle->SetOptStat(0);

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

void plotGraph(std::string name, std::string xTitle,
		std::vector<double> xValues, double xError, std::vector<double> yValues,
		std::vector<double> yErrors, std::string subdir, double fitRangeStart,
		double fitRangeEnd) {

	double xErrors[xValues.size()];
	for (unsigned int i = 0; i < xValues.size(); i++) {
		xErrors[i] = xError;
	}

	TGraphErrors hitWidthVsV(xValues.size(), &xValues[0], &yValues[0], xErrors,
			&yErrors[0]);
	hitWidthVsV.SetTitle(name.c_str());
	hitWidthVsV.SetName(hitWidthVsV.GetTitle());
	hitWidthVsV.GetXaxis()->SetTitle(xTitle.c_str());
	hitWidthVsV.GetYaxis()->SetTitle("average hit width [strips]");
	hitWidthVsV.SetDrawOption("AP");
	hitWidthVsV.SetMarkerColor(4);
	hitWidthVsV.SetMarkerStyle(21);

	std::cout << "########################################" << std::endl;
	std::cout << "Fitting " << name << std::endl;
	std::cout << "########################################" << std::endl;
	TF1 f1("f1", "pol1", fitRangeStart, fitRangeEnd);
	//hitWidthVsV.Fit("pol1", "", "", fitRangeStart, fitRangeEnd);
	hitWidthVsV.Fit(&f1, "R");
	hitWidthVsV.Write(hitWidthVsV.GetTitle());

	writeToPdf<TGraph>(&hitWidthVsV, subdir, "AP");
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
			yValueErrorsFiltered, subfolder.str(), fitRangeStart, fitRangeEnd);
}
