#ifndef MMQuickEvent_H 
#define MMQuickEvent_H

#include "CCommonIncludes.h"
#include "MapFile.h"

using namespace std;

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

template<class T1, class T2>
struct sort_pair_first {
	bool operator()(const std::pair<T1, T2>&left,
			const std::pair<T1, T2>&right) {
		return left.first < right.first;
	}
};

struct cutStatistics_t {
	TH1F* timingCuts;
	TH1F* chargeCuts;
	TH1F* timeCoincidenceCuts;
	TH1F* absolutePositionCuts;
	TH1F* proportionXCuts;
	TH1F* proportionYCuts;
	TH1F* fitMeanMaxChargeDistanceCuts;
	TH1F* fitProblemCuts;
};

class MMQuickEvent {
public:
	MMQuickEvent(vector<string> vecFilenames, string Tree_Name,
			int NumberOfEvents = 10) {
		m_tchain = new TChain(Tree_Name.c_str());
		for (unsigned int i = 0; i < vecFilenames.size(); i++) {
			m_tchain->Add(vecFilenames[i].c_str());
			m_tchain->AddFriend("data", vecFilenames[i].c_str());
		}
		cleanVariables();
		addBranches();
		m_actEventNumber = 0;
		m_NumberOfEvents = m_tchain->GetEntries();
		if (NumberOfEvents != -1)
			m_NumberOfEvents = NumberOfEvents;
		cout << "[MMQuickEvent] Number of Events loaded: " << m_NumberOfEvents
				<< endl;

		maxChargeX = 0;
		stripWithMaxChargeX = 0;
		timeSliceOfMaxChargeX = 0;
		maxChargeY = 0;
		stripWithMaxChargeY = 0;
		timeSliceOfMaxChargeY = 0;

		numberOfXHits = 0;
		numberOfYHits = 0;
	}

	bool getNextEvent() {
		if (m_actEventNumber >= m_NumberOfEvents) {
			cout << endl;
			return false;
		}
		if (m_actEventNumber == 0)
			cout << "[MMQuickEvent] Looping over Events" << endl;
		if (m_actEventNumber % (m_NumberOfEvents / 100) == 0) {
			cout << '\r' << "[MMQuickEvent] "
					<< TMath::Nint(
							m_actEventNumber / ((float) m_NumberOfEvents)
									* 100.) << "% done (Event "
					<< m_actEventNumber << "/" << m_NumberOfEvents << ")...";
			cout.flush();
		} else if (m_actEventNumber == m_NumberOfEvents - 1) {
			cout << '\r'
					<< "[MMQuickEvent] Done!                              ";
			cout.flush();
		}
		m_tchain->GetEvent(m_actEventNumber);
		m_actEventNumber++;

		maxChargeX = 0;
		stripWithMaxChargeX = 0;
		timeSliceOfMaxChargeX = 0;
		maxChargeY = 0;
		stripWithMaxChargeY = 0;
		timeSliceOfMaxChargeY = 0;

		numberOfXHits = 0;
		numberOfYHits = 0;

		return true;
	}

	void cleanVariables() {
		apv_fecNo = 0;
		apv_id = 0;
		apv_ch = 0;
		mm_id = 0;
		mm_readout = 0;
		mm_strip = 0;
		apv_q = 0;

		apv_qmax = 0;
		apv_tbqmax = 0;
	}

	int getEventNumber() {
		return m_NumberOfEvents;
	}

	int getCurrentEventNumber() {
		return m_actEventNumber;
	}

	void addBranches() {
		m_tchain->SetBranchAddress("apv_evt", &apv_evt);
		m_tchain->SetBranchAddress("time_s", &time_s);
		m_tchain->SetBranchAddress("time_us", &time_us);
		m_tchain->SetBranchAddress("apv_fecNo", &apv_fecNo);
		m_tchain->SetBranchAddress("apv_id", &apv_id);
		m_tchain->SetBranchAddress("apv_ch", &apv_ch);
		m_tchain->SetBranchAddress("mm_id", &mm_id);
		m_tchain->SetBranchAddress("mm_readout", &mm_readout);
		m_tchain->SetBranchAddress("mm_strip", &mm_strip);
		m_tchain->SetBranchAddress("apv_q", &apv_q);
		m_tchain->SetBranchAddress("apv_presamples", &apv_presamples);

		m_tchain->SetBranchAddress("apv_qmax", &apv_qmax);
		m_tchain->SetBranchAddress("apv_tbqmax", &apv_tbqmax);
	}

	// functions to select if hit is in X or Y according to APV ID and mapping while data acquisition
	static bool isX(int id) {
		if (id == APVIDMM_X0 || id == APVIDMM_X1 || id == APVIDMM_X2) {
			return true;
		} else {
			return false;
		}
	}

	static bool isY(int id) {
		if (id == APVIDMM_Y0 || id == APVIDMM_Y1 || id == APVIDMM_Y2) {
			return true;
		} else {
			return false;
		}
	}

public:
	TChain *m_tchain;
	int m_actEventNumber;
	int m_NumberOfEvents;

	/// Event Information
	// Declaration of leaf types
	UInt_t apv_evt;
	Int_t time_s;
	Int_t time_us;
	vector<unsigned int> *apv_fecNo;
	vector<unsigned int> *apv_id;
	vector<unsigned int> *apv_ch;
	vector<string> *mm_id;
	vector<unsigned int> *mm_readout;
	vector<unsigned int> *mm_strip;
	vector<vector<short> > *apv_q;
	UInt_t apv_presamples;

	vector<short> *apv_qmax;
	vector<short> *apv_tbqmax;

	short maxChargeX;
	int stripWithMaxChargeX;
	int timeSliceOfMaxChargeX;
	short maxChargeY;
	int stripWithMaxChargeY;
	int timeSliceOfMaxChargeY;

	unsigned short numberOfXHits;
	unsigned short numberOfYHits;

	vector<std::pair<unsigned int, short> > stripAndChargeAtMaxChargeTimeX; // absolute strip number and charges of all strips at fixed time slice (being the maximum charge time)
	vector<std::pair<unsigned int, short> > stripAndChargeAtMaxChargeTimeY;

	/**
	 * Returns true if the neighbour strips of the strip with maximum charge are within a given range
	 */
	void generateFixedTimeCrossSections() {
		/*
		 * Store the charge values of every strip number for the time slice with
		 * the maximum charge found in one event for X and Y separately (cross section
		 * for time slices with max charge)
		 */
		stripAndChargeAtMaxChargeTimeX.clear();
		stripAndChargeAtMaxChargeTimeY.clear();

		// Iterate through all strips
		for (unsigned int strip = 0; strip != (*apv_q).size(); strip++) {
			unsigned int apvID = (*apv_id)[strip];
			if (MMQuickEvent::isX(apvID)) { // X axis
				stripAndChargeAtMaxChargeTimeX.push_back(
						std::make_pair((*mm_strip)[strip]/*Strip number*/,
								(*apv_q)[strip][timeSliceOfMaxChargeX]/*Charge*/));
			} else { // Y axis
				stripAndChargeAtMaxChargeTimeY.push_back(
						std::make_pair((*mm_strip)[strip]/*Strip number*/,
								(*apv_q)[strip][timeSliceOfMaxChargeY]/*Charge*/));
			}
		}
	}

	bool runProportionCut(TH2F* maxNeighbourHisto,
			vector<std::pair<unsigned int, short> > stripAndChargeAtMaxChargeTime,
			short maxCharge, std::vector<std::pair<int, int> > proportionLimits,
			cutStatistics_t& cutStat, TH1F* proportionCuts) {

		// Sort cross section by absolute strip numbers (first entry in pairs)
		std::sort(stripAndChargeAtMaxChargeTime.begin(),
				stripAndChargeAtMaxChargeTime.end(),
				sort_pair_first<unsigned int, short>());

		/*
		 * Now the array is sorted, the position of the maximal charge strip is unknown -> search for it again
		 */
		unsigned int positionOfMaxCharge;

		for (positionOfMaxCharge = 0;
				positionOfMaxCharge < stripAndChargeAtMaxChargeTime.size();
				positionOfMaxCharge++) {
			if (stripAndChargeAtMaxChargeTime[positionOfMaxCharge].second
					== maxCharge) {
				break;
			}
		}

		/*
		 * Start at strip N left to the maximum charge strip where N is the number of entries in the propotion limits file
		 */
		bool accepEvent = true;
		const int maxDistance = proportionLimits.size();
		for (int deltaStrip = -maxDistance; deltaStrip <= maxDistance;
				deltaStrip++) {
			if (deltaStrip == 0) {
				// don't store the maximum strip itself
				continue;
			}

			// Look at the x/y strips deltaStrip away from the maximal charge strip in x/y
			int stripIndexInCrossSection = positionOfMaxCharge + deltaStrip;

			// Check if all neighbour strips are available
			bool tooFarToTheLeft =
					stripAndChargeAtMaxChargeTime[positionOfMaxCharge].first
							+ deltaStrip <= 0; // absolute strip too far to the left
			bool tooFarToTheRight =
					stripAndChargeAtMaxChargeTime[positionOfMaxCharge].first
							+ deltaStrip > xStrips; // absolute strip too far to the left

			bool lowerLimitIsZero = proportionLimits[abs(deltaStrip) - 1].first
					== 0; // only check if strip has charge stored if lower limit is not zero

			bool stripChargeIsStored =
					(stripIndexInCrossSection >= 0 // too far to the left in the array?
							&& stripIndexInCrossSection
									< stripAndChargeAtMaxChargeTime.size() // too far to the right in the array
							&& stripAndChargeAtMaxChargeTime[stripIndexInCrossSection].first // check if there is a jump in the array
									== stripAndChargeAtMaxChargeTime[positionOfMaxCharge].first
											+ deltaStrip);

			double proportion = NAN;
			if (tooFarToTheRight || tooFarToTheLeft) {
				if (accepEvent) {
					cutStat.absolutePositionCuts->Fill(1);
				} else {
					cutStat.absolutePositionCuts->Fill(0);
				}
				accepEvent = false;
			} else if (!stripChargeIsStored && !lowerLimitIsZero) {
				if (accepEvent) {
					proportionCuts->Fill(1);
				} else {
					proportionCuts->Fill(0);
				}
				accepEvent = false;
			} else {
				if (stripChargeIsStored) {
					proportion =
							100
									* stripAndChargeAtMaxChargeTime[stripIndexInCrossSection].second
									/ (double) maxCharge;

					// Cut for neighbours of maximum bin
					if (proportion < proportionLimits[abs(deltaStrip) - 1].first
							|| proportion
									> proportionLimits[abs(deltaStrip) - 1].second) {
						if (accepEvent) {
							proportionCuts->Fill(1);
						} else {
							proportionCuts->Fill(0);
						}
						accepEvent = false;
					}
				}
			}

			if (proportion != NAN) {
				// Fill histogramm mmhitneighboursX and mmhitneighboursY
				maxNeighbourHisto->Fill((deltaStrip), proportion);
			}
		}
		return accepEvent;
	}

	void findMaxCharge() {

		vector<unsigned int> apvIDofStrip = *apv_id; // isX(apvIDofStrip[i]) returns true if the i-th strip is X-layer
		vector<unsigned int> stripNumShowingSignal = *mm_strip; // stripNumShowingSignal[i] is absolute strip number (strips without charge are not stored anywhere)
		vector<vector<short> > chargeOfStripOfTime = *apv_q; // chargeOfStripOfTime[i][j] is the charge of strip i in time slice j (matrix of whole event)
		vector<short> maxChargeOfStrip = *apv_qmax; // maxChargeOfStrip[i] is the maxmimal measured charge of strip i of all time slices
		vector<short> timeSliceOfMaxChargeOfStrip = *apv_tbqmax; // timeSliceOfMaxChargeOfStrip[i] is the time slice of the corresponding maximum charge (see above)

		/*
		 * Iterate through all strips and check if it is X or Y data. Compare the maximum charge
		 * of the strip with the maximum charge found so far for the current axis. Store current charge, strip number
		 * and time slice with the maximum charge if the current charge is larger than before.
		 */
		for (unsigned int strip = 1; strip != maxChargeOfStrip.size();
				strip++) {
			unsigned int apvID = apvIDofStrip[strip];
			if (isX(apvID)) { // X axis
				numberOfXHits++;
				if (maxChargeOfStrip[strip] > maxChargeX) {
					maxChargeX = maxChargeOfStrip[strip];
					stripWithMaxChargeX = strip;
					timeSliceOfMaxChargeX = timeSliceOfMaxChargeOfStrip[strip];
				}
			} else { // Y axis
				numberOfYHits++;
				if (maxChargeOfStrip[strip] > maxChargeY) {
					maxChargeY = maxChargeOfStrip[strip];
					stripWithMaxChargeY = strip;
					timeSliceOfMaxChargeY = timeSliceOfMaxChargeOfStrip[strip];
				}
			}
		}
	}
};

#endif
