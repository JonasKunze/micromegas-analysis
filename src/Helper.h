/*
 * Helper.cxx
 *
 *  Created on: Feb 8, 2015
 *      Author: kunzejo
 */
#ifndef HELPER_H
#define HELPER_H

#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <cstdlib>
#include <sstream>
#include <string>

class TH2F;

//set output path and name of output files
const std::string inPath = "/localscratch/praktikum/data/";	//Path of the Input
const std::string outPath = "/localscratch/praktikum/output/"; // Path of the Output
//const string outPath = "/tmp/output/"; // Path of the Output
const std::string appendName = "";				// Name of single measurements
const std::string combinedPlotsFile = "combined.root";// Name of the file for the combined results of all runs (hier muss jeder Tag einzeln analysiert werden! Da Zeile 79-84(driftStart...ampSteps) für jeden Tag anders war.

template<typename T>
void writeToPdf(T* object, std::string subfolder, std::string drawOptions) {
	gStyle->SetOptFit(1111);

	std::string histoName(object->GetName());
	gStyle->SetStatW(0.2);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.47);


	if(dynamic_cast<TH1*>(object)){
		gStyle->SetStatY(0.9);
		gStyle->SetStatX(0.9);
	}
	if (histoName == "timeDistribution") {
		gStyle->SetStatY(0.5);
		gStyle->SetStatX(0.6);
	}
	std::stringstream pdfName;
	pdfName << outPath << "/" << subfolder << "/";

	std::stringstream mkdir;
	mkdir << "mkdir -p " << pdfName.str();
	system(mkdir.str().c_str());

	pdfName << object->GetName() << ".pdf";

	TCanvas canvas("c", "data", 200, 10, 700, 500);
	canvas.SetLeftMargin(0.11);
	object->GetYaxis()->SetTitleOffset(1.5);
	object->Draw(drawOptions.c_str());
	canvas.Print(pdfName.str().c_str(), "pdf");
}
void writeTH2FToPdf(TH2F* object, std::string subfolder,
		std::string drawOptions);
#endif
