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
		m_mapFile["VD50VA500"] = new TFile(
				(path + appendName + "_VD50VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD125VA500"] = new TFile(
				(path + appendName + "_VD125VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD200VA500"] = new TFile(
				(path + appendName + "_VD200VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD275VA500"] = new TFile(
				(path + appendName + "_VD275VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD350VA500"] = new TFile(
				(path + appendName + "_VD350VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD50VA525"] = new TFile(
				(path + appendName + "_VD50VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD125VA525"] = new TFile(
				(path + appendName + "_VD125VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD200VA525"] = new TFile(
				(path + appendName + "_VD200VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD275VA525"] = new TFile(
				(path + appendName + "_VD275VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD350VA525"] = new TFile(
				(path + appendName + "_VD350VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD50VA550"] = new TFile(
				(path + appendName + "_VD50VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD125VA550"] = new TFile(
				(path + appendName + "_VD125VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD200VA550"] = new TFile(
				(path + appendName + "_VD200VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD275VA550"] = new TFile(
				(path + appendName + "_VD275VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD350VA550"] = new TFile(
				(path + appendName + "_VD350VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD172VA500"] = new TFile(
				(path + appendName + "_VD172VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD430VA500"] = new TFile(
				(path + appendName + "_VD430VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD688VA500"] = new TFile(
				(path + appendName + "_VD688VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD947VA500"] = new TFile(
				(path + appendName + "_VD947VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD1205VA500"] = new TFile(
				(path + appendName + "_VD1205VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD172VA525"] = new TFile(
				(path + appendName + "_VD172VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD430VA525"] = new TFile(
				(path + appendName + "_VD430VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD688VA525"] = new TFile(
				(path + appendName + "_VD688VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD947VA525"] = new TFile(
				(path + appendName + "_VD947VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD1205VA525"] = new TFile(
				(path + appendName + "_VD1205VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD172VA550"] = new TFile(
				(path + appendName + "_VD172VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD430VA550"] = new TFile(
				(path + appendName + "_VD430VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD688VA550"] = new TFile(
				(path + appendName + "_VD688VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD947VA550"] = new TFile(
				(path + appendName + "_VD947VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD117VA500"] = new TFile(
				(path + appendName + "_VD117VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD291VA500"] = new TFile(
				(path + appendName + "_VD291VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD467VA500"] = new TFile(
				(path + appendName + "_VD467VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD642VA500"] = new TFile(
				(path + appendName + "_VD642VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD817VA500"] = new TFile(
				(path + appendName + "_VD817VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD117VA525"] = new TFile(
				(path + appendName + "_VD117VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD291VA525"] = new TFile(
				(path + appendName + "_VD291VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD467VA525"] = new TFile(
				(path + appendName + "_VD467VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD642VA525"] = new TFile(
				(path + appendName + "_VD642VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD817VA525"] = new TFile(
				(path + appendName + "_VD817VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD117VA550"] = new TFile(
				(path + appendName + "_VD117VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD291VA550"] = new TFile(
				(path + appendName + "_VD291VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD467VA550"] = new TFile(
				(path + appendName + "_VD467VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD642VA550"] = new TFile(
				(path + appendName + "_VD642VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD817VA550"] = new TFile(
				(path + appendName + "_VD817VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD89VA500"] = new TFile(
				(path + appendName + "_VD89VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD222VA500"] = new TFile(
				(path + appendName + "_VD222VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD355VA500"] = new TFile(
				(path + appendName + "_VD355VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD488VA500"] = new TFile(
				(path + appendName + "_VD488VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD622VA500"] = new TFile(
				(path + appendName + "_VD622VA500.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD89VA525"] = new TFile(
				(path + appendName + "_VD89VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD222VA525"] = new TFile(
				(path + appendName + "_VD222VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD355VA525"] = new TFile(
				(path + appendName + "_VD355VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD488VA525"] = new TFile(
				(path + appendName + "_VD488VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD622VA525"] = new TFile(
				(path + appendName + "_VD622VA525.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD89VA550"] = new TFile(
				(path + appendName + "_VD89VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD222VA550"] = new TFile(
				(path + appendName + "_VD222VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD355VA550"] = new TFile(
				(path + appendName + "_VD355VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD488VA550"] = new TFile(
				(path + appendName + "_VD488VA550.root").c_str(),
				(Option_t*) "RECREATE");
		m_mapFile["VD622VA550"] = new TFile(
				(path + appendName + "_VD622VA550.root").c_str(),
				(Option_t*) "RECREATE");
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

	vector<string> getFileName(string type) {
		vector<string> vec_filename;
		cout
				<< "======================================================================"
				<< endl;
		cout << "================== Processing dataset of : " << type
				<< " ====================" << endl;
		cout
				<< "======================================================================"
				<< endl;

		if (type == "VD50VA500") {
			vec_filename.push_back(data_dir + "run737.root");
		} else if (type == "VD125VA500") {
			vec_filename.push_back(data_dir + "run376.root");
		} else if (type == "VD200VA500") {
			vec_filename.push_back(data_dir + "run378.root");
		} else if (type == "VD275VA500") {
			vec_filename.push_back(data_dir + "run380.root");
		} else if (type == "VD350VA500") {
			vec_filename.push_back(data_dir + "run382.root");
		} else if (type == "VD50VA525") {
			vec_filename.push_back(data_dir + "run384.root");
		} else if (type == "VD125VA525") {
			vec_filename.push_back(data_dir + "run387.root");
		} else if (type == "VD200VA525") {
			vec_filename.push_back(data_dir + "run389.root");
		} else if (type == "VD275VA525") {
			vec_filename.push_back(data_dir + "run391.root");
		} else if (type == "VD350VA525") {
			vec_filename.push_back(data_dir + "run393.root");
		} else if (type == "VD50VA550") {
			vec_filename.push_back(data_dir + "run395.root");
		} else if (type == "VD125VA550") {
			vec_filename.push_back(data_dir + "run397.root");
		} else if (type == "VD200VA550") {
			vec_filename.push_back(data_dir + "run399.root");
		} else if (type == "VD275VA550") {
			vec_filename.push_back(data_dir + "run401.root");
		} else if (type == "VD350VA550") {
			vec_filename.push_back(data_dir + "run407.root");
		} else if (type == "VD172VA500") {
			vec_filename.push_back(data_dir + "run411.root");
		} else if (type == "VD430VA500") {
			vec_filename.push_back(data_dir + "run413.root");
		} else if (type == "VD688VA500") {
			vec_filename.push_back(data_dir + "run415.root");
		} else if (type == "VD947VA500") {
			vec_filename.push_back(data_dir + "run417.root");
		} else if (type == "VD1205VA500") {
			vec_filename.push_back(data_dir + "run419.root");
		} else if (type == "VD172VA525") {
			vec_filename.push_back(data_dir + "run424.root");
		} else if (type == "VD430VA525") {
			vec_filename.push_back(data_dir + "run426.root");
		} else if (type == "VD688VA525") {
			vec_filename.push_back(data_dir + "run428.root");
		} else if (type == "VD947VA525") {
			vec_filename.push_back(data_dir + "run430.root");
		} else if (type == "VD1205VA525") {
			vec_filename.push_back(data_dir + "run432.root");
		} else if (type == "VD172VA550") {
			vec_filename.push_back(data_dir + "run434.root");
		} else if (type == "VD430VA550") {
			vec_filename.push_back(data_dir + "run442.root");
		} else if (type == "VD688VA550") {
			vec_filename.push_back(data_dir + "run445.root");
		} else if (type == "VD947VA550") {
			vec_filename.push_back(data_dir + "run450.root");
		} else if (type == "VD117VA500") {
			vec_filename.push_back(data_dir + "run453.root");
		} else if (type == "VD291VA500") {
			vec_filename.push_back(data_dir + "run455.root");
		} else if (type == "VD467VA500") {
			vec_filename.push_back(data_dir + "run457.root");
		} else if (type == "VD642VA500") {
			vec_filename.push_back(data_dir + "run459.root");
		} else if (type == "VD817VA500") {
			vec_filename.push_back(data_dir + "run461.root");
		} else if (type == "VD117VA525") {
			vec_filename.push_back(data_dir + "run463.root");
		} else if (type == "VD291VA525") {
			vec_filename.push_back(data_dir + "run465.root");
		} else if (type == "VD467VA525") {
			vec_filename.push_back(data_dir + "run467.root");
		} else if (type == "VD642VA525") {
			vec_filename.push_back(data_dir + "run469.root");
		} else if (type == "VD817VA525") {
			vec_filename.push_back(data_dir + "run472.root");
		} else if (type == "VD117VA550") {
			vec_filename.push_back(data_dir + "run475.root");
		} else if (type == "VD291VA550") {
			vec_filename.push_back(data_dir + "run477.root");
		} else if (type == "VD467VA550") {
			vec_filename.push_back(data_dir + "run480.root");
		} else if (type == "VD642VA550") {
			vec_filename.push_back(data_dir + "run482.root");
		} else if (type == "VD817VA550") {
			vec_filename.push_back(data_dir + "run484.root");
		} else if (type == "VD89VA500") {
			vec_filename.push_back(data_dir + "run487.root");
		} else if (type == "VD222VA500") {
			vec_filename.push_back(data_dir + "run490.root");
		} else if (type == "VD355VA500") {
			vec_filename.push_back(data_dir + "run492.root");
		} else if (type == "VD488VA500") {
			vec_filename.push_back(data_dir + "run494.root");
		} else if (type == "VD622VA500") {
			vec_filename.push_back(data_dir + "run496.root");
		} else if (type == "VD89VA525") {
			vec_filename.push_back(data_dir + "run505.root");
		} else if (type == "VD222VA525") {
			vec_filename.push_back(data_dir + "run507.root");
		} else if (type == "VD355VA525") {
			vec_filename.push_back(data_dir + "run509.root");
		} else if (type == "VD488VA525") {
			vec_filename.push_back(data_dir + "run511.root");
		} else if (type == "VD622VA525") {
			vec_filename.push_back(data_dir + "run513.root");
		} else if (type == "VD89VA550") {
			vec_filename.push_back(data_dir + "run515.root");
		} else if (type == "VD222VA550") {
			vec_filename.push_back(data_dir + "run518.root");
		} else if (type == "VD355VA550") {
			vec_filename.push_back(data_dir + "run520.root");
		} else if (type == "VD488VA550") {
			vec_filename.push_back(data_dir + "run522.root");
		} else if (type == "VD622VA550") {
			vec_filename.push_back(data_dir + "run524.root");
		}

		return vec_filename;
	}

private:
	map<string, TFile*> m_mapFile;
	string data_dir; // was "../../PhD/Detector/micromega_data/" before
	string path;
	string appendName;
};

#endif
