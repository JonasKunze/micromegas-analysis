/*
 * Helper.cxx
 *
 *  Created on: Feb 8, 2015
 *      Author: kunzejo
 */

#include "Helper.h"

void writeTH2FToPdf(TH2F* object, std::string subfolder,
		std::string drawOptions) {
	std::stringstream pdfName;
	pdfName << outPath << "/" << subfolder << "/";

	std::stringstream mkdir;
	mkdir << "mkdir -p " << pdfName.str();
	system(mkdir.str().c_str());

	pdfName << object->GetName() << ".pdf";

	TCanvas canvas("c", "data", 200, 10, 700, 500);
	canvas.SetRightMargin(0.15);
	object->GetZaxis()->SetTitleOffset(1.2);
	object->Draw(drawOptions.c_str());
	canvas.Print(pdfName.str().c_str(), "pdf");
}

