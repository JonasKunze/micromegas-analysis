/*
 * CutStatistics.cxx
 *
 *  Created on: Feb 7, 2015
 *      Author: kunzejo
 */

#include "CutStatistic.h"

#include "MMQuickEvent.h"

#define MAX_EVENT_DISPLAYS_PER_CUT 5

std::vector<CutStatistic*> CutStatistic::instances;

void CutStatistic::Fill(double value, MMQuickEvent* event, std::string suffix) {
	counterHistogram.Fill(value);

	if (value == 1
			&& eventDisplaysCut.size() < MAX_EVENT_DISPLAYS_PER_CUT * 2) {
		TH2F *eventDisplayX, *eventDisplayY;

		std::stringstream eventDisplaySuffix;
		eventDisplaySuffix << "-" << counterHistogram.GetName() << suffix;
		event->generateEventDisplay(eventDisplayX, eventDisplayY,
				eventDisplaySuffix.str());
		eventDisplaysCut.push_back(eventDisplayX);
		eventDisplaysCut.push_back(eventDisplayY);
	} else if (value == 0
			&& eventDisplaysAccepted.size() < MAX_EVENT_DISPLAYS_PER_CUT * 2) {
		TH2F *eventDisplayX, *eventDisplayY;

		std::stringstream eventDisplaySuffix;
		eventDisplaySuffix << "-" << counterHistogram.GetName() << suffix;
		event->generateEventDisplay(eventDisplayX, eventDisplayY,
				eventDisplaySuffix.str());
		eventDisplaysAccepted.push_back(eventDisplayX);
		eventDisplaysAccepted.push_back(eventDisplayY);
	}
}

