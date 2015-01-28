#!/bin/sh
#cat runs.txt | grep phy | awk '{print "m_mapFile[\"VD"$3"VA"$2"\"] = new TFile((path+appendName+\"_VD"$3"VA"$2".root\").c_str(), (Option_t*)\"RECREATE\");"}'
cat runs.txt | grep phy | awk '{print "if(type==\"VD"$3"VA"$2"\"){vec_filename.push_back(data_dir + \"run"$1".root\");} else"}'
