#ifndef MMQuickEvent_H 
#define MMQuickEvent_H

#include "CCommonIncludes.h"
using namespace std;

class MMQuickEvent {
	public:
		MMQuickEvent(vector<string> vecFilenames, string Tree_Name, int NumberOfEvents=10) {
			m_tchain = new TChain(Tree_Name.c_str());
			for (unsigned int i=0; i<vecFilenames.size(); i++) {
				m_tchain->Add(vecFilenames[i].c_str());
				m_tchain->AddFriend("data", vecFilenames[i].c_str());
			}
			cleanVariables();
			addBranches();
			m_actEventNumber = 0;
			m_NumberOfEvents = m_tchain->GetEntries();
			if (NumberOfEvents!=-1) m_NumberOfEvents = NumberOfEvents;
			cout<<"[MMQuickEvent] Number of Events loaded: "<<m_NumberOfEvents<<endl;
		}

		bool getNextEvent() {
			if (m_actEventNumber >= m_NumberOfEvents) {
				cout << endl;
				return false;
			}
			if (m_actEventNumber == 0) cout << "[MMQuickEvent] Looping over Events" << endl;
			if (m_actEventNumber%(m_NumberOfEvents/100) == 0) {
				cout << '\r' << "[MMQuickEvent] " << TMath::Nint(m_actEventNumber/((float)m_NumberOfEvents)*100.) << "% done (Event " << m_actEventNumber << "/" << m_NumberOfEvents << ")...";
				cout.flush();
			} else if (m_actEventNumber == m_NumberOfEvents-1) {
				cout << '\r' << "[MMQuickEvent] Done!                              ";
				cout.flush();
			}
			m_tchain->GetEvent(m_actEventNumber);
			m_actEventNumber++;
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

		double getEventNumber() {
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

	public:
		TChain 					*m_tchain;
		int 					m_actEventNumber;
		int 					m_NumberOfEvents;

		/// Event Information
		// Declaration of leaf types
		UInt_t					apv_evt;
		Int_t					time_s;
		Int_t					time_us;
		vector<unsigned int> 	*apv_fecNo;
		vector<unsigned int> 	*apv_id;
		vector<unsigned int> 	*apv_ch;
		vector<string> 			*mm_id;
		vector<unsigned int> 	*mm_readout;
		vector<unsigned int> 	*mm_strip;
		vector<vector<short> > 	*apv_q;
		UInt_t 					apv_presamples;

		vector<short>			*apv_qmax;
		vector<short>			*apv_tbqmax;
};

#endif