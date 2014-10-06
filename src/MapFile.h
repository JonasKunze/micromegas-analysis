#ifndef MapFile_H 
#define MapFile_H

#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
using namespace std;

class MapFile {
private:
	void createFile() {
		m_mapFile["VD100VA500"] = new TFile((path+appendName+"_VD100VA500.root").c_str(), (Option_t*)"RECREATE");		// string has to be like VDvalueVAvalue, value must have 3 digits
		m_mapFile["VD100VA525"] = new TFile((path+appendName+"_VD100VA500.root").c_str(), (Option_t*)"RECREATE");
	}

public:
		/**
		 * Consctructor
		 */
		MapFile(string data_dir, string path, string appendName) {
			this->data_dir = data_dir;
			this->path = path;
			this->appendName = appendName;
			createFile();
		}

		map<string, TFile*> getFile() {
			return m_mapFile;
		}

		vector<string> getFileName (string type) {
			vector<string> vec_filename;
			cout<<"======================================================================"<<endl;
			cout<<"================== Processing dataset of : " << type << " ===================="<<endl;
			cout<<"======================================================================"<<endl;

			if (type == "VD100VA500") {
				vec_filename.push_back(data_dir + "run332.root");
			} else if (type == "VD100VA525") {
				vec_filename.push_back(data_dir + "run334.root");

			}

			return vec_filename;

		}


	private:
		map<string, TFile*>	m_mapFile;
		string 			data_dir; // was "../../PhD/Detector/micromega_data/" before
		string 			path;
		string			appendName;
};

#endif
