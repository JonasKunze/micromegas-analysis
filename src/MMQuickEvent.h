#ifndef MMQuickEvent_H 
#define MMQuickEvent_H

#include "CCommonIncludes.h"

using namespace std;

// Mapping of APV-Chips
const int APVIDMM_X0 = 5;
const int APVIDMM_X1 = 4;
const int APVIDMM_X2 = 6;
const int APVIDMM_Y0 = 0;
const int APVIDMM_Y1 = 1;
const int APVIDMM_Y2 = 2;

const int NumberOfMaxHitNeighboursToBeStored = 5;

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

	vector<short> chargeOfStripAtMaxChargeTimeX; // charges of all strips at fixed time slice (being the maximum charge time)
	vector<short> chargeOfStripAtMaxChargeTimeY;
	vector<unsigned int> numberOfStripAtMaxChargeTimeX; // absolute strip number of all strips
	vector<unsigned int> numberOfStripAtMaxChargeTimeY;

	void generateFixTimeCrossSection() {
		/*
		 * Store the charge values of every strip number for the time slice with
		 * the maximum charge found in one event for X and Y separately (cross section
		 * for time slices with max charge)
		 */

		chargeOfStripAtMaxChargeTimeX.clear();
		chargeOfStripAtMaxChargeTimeY.clear();
		numberOfStripAtMaxChargeTimeX.clear();
		numberOfStripAtMaxChargeTimeY.clear();

		// Iterate through all strips
		for (unsigned int strip = 0; strip != (*apv_q).size(); strip++) {
			unsigned int apvID = (*apv_id)[strip];
			if (MMQuickEvent::isX(apvID)) { // X axis
				chargeOfStripAtMaxChargeTimeX.push_back(
						(*apv_q)[strip][timeSliceOfMaxChargeX]);
				numberOfStripAtMaxChargeTimeX.push_back((*mm_strip)[strip]);
			} else { // Y axis
				chargeOfStripAtMaxChargeTimeY.push_back(
						(*apv_q)[strip][timeSliceOfMaxChargeY]);
				numberOfStripAtMaxChargeTimeY.push_back((*mm_strip)[strip]);
			}
		}
	}

	void findMaxCharge(TH2F* maxNeighbourHisto) {

		vector<unsigned int> apvIDofStrip = *apv_id; // isX(apvIDofStrip[i]) returns true if the i-th strip is X-layer
		vector<unsigned int> stripNumShowingSignal = *mm_strip; // stripNumShowingSignal[i] is absolute strip number (strips without charge are not stored anywhere)
		vector<vector<short> > chargeOfStripOfTime = *apv_q; // chargeOfStripOfTime[i][j] is the charge of strip i in time slice j (matrix of whole event)
		vector<short> maxChargeOfStrip = *apv_qmax; // maxChargeOfStrip[i] is the maxmimal measured charge of strip i of all time slices
		vector<short> timeSliceOfMaxChargeOfStrip = *apv_tbqmax; // timeSliceOfMaxChargeOfStrip[i] is the time slice of the corresponding maximum charge (see above)

		std::map<unsigned int, unsigned short> chargesOfAllStripsX;
		std::map<unsigned int, unsigned short> chargesOfAllStripsY;
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
				chargesOfAllStripsX[stripNumShowingSignal[strip]] =
						maxChargeOfStrip[strip];
			} else { // Y axis
				numberOfYHits++;
				if (maxChargeOfStrip[strip] > maxChargeY) {
					maxChargeY = maxChargeOfStrip[strip];
					stripWithMaxChargeY = strip;
					timeSliceOfMaxChargeY = timeSliceOfMaxChargeOfStrip[strip];
				}
				chargesOfAllStripsY[stripNumShowingSignal[strip]] =
						maxChargeOfStrip[strip];
			}
		}

		for (int delta = -NumberOfMaxHitNeighboursToBeStored;
				delta < NumberOfMaxHitNeighboursToBeStored; delta++) {

			if (delta == 0) {
				// don't store the maxim strip itself
				continue;
			}

			int stripX = timeSliceOfMaxChargeOfStrip[stripWithMaxChargeX]
					+ delta;
			int stripY = timeSliceOfMaxChargeOfStrip[stripWithMaxChargeY]
					+ delta;

			if (chargesOfAllStripsX.find(stripX) != chargesOfAllStripsX.end()) {
				maxNeighbourHisto->Fill(abs(delta),
						100*(double)chargesOfAllStripsX[stripX]/maxChargeX );

			}
			if (chargesOfAllStripsY.find(stripY) != chargesOfAllStripsY.end()) {
				maxNeighbourHisto->Fill(abs(delta),
						100*(double)chargesOfAllStripsY[stripY]/maxChargeY);
			}
		}
	}
};

#endif
