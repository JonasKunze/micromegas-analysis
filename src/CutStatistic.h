/*
 * CutStatistics.h
 *
 *  Created on: Feb 7, 2015
 *      Author: kunzejo
 */

#ifndef CUTSTATISTIC_H_
#define CUTSTATISTIC_H_

#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <vector>

class MMQuickEvent;

class CutStatistic {
public:
	static std::vector<CutStatistic*> instances;
	TH1F counterHistogram;
	std::vector<TH2F*> eventDisplaysCut;
	std::vector<TH2F*> eventDisplaysAccepted;

	CutStatistic(std::string name):counterHistogram(name.c_str(), ";accepted/cut ;entries", 2,
			-0.5, 1.5) {
		instances.push_back(this);
	}

	void Fill(double value, MMQuickEvent* event);

	const char* getName() {
		return counterHistogram.GetName();
	}
};

#endif /* CUTSTATISTIC_H_ */
