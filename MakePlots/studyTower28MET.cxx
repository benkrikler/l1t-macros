#include "../Plotting/TL1RootHists.h"
#include "../Plotting/TL1Rates.h"
#include "../Plotting/TL1Turnon.h"
#include "../Event/TL1EventClass.h"
#include "../Config/ntuple_cfg.h"

#include <string>
#include <iostream>
#include <TH1F.h>

using std::cout;
using std::endl;

typedef std::map<std::string,TL1Plots*> PlotList;
PlotList MakePlots(ntuple_cfg*);
void FillPlots(int i_entry, PlotList& plots, const TL1EventClass*);

PlotList MakePlots(ntuple_cfg* dataset){
    PlotList plots;

    std::string outName = dataset->triggerName;

    double met_max=100;
    int met_bins=50;
    int phi_bins=72;

    // Tower 28 MET
    std::string name = "tower28MET";
    TH1* tower28MET_template=new TH1F(name.c_str(),"#splitline{MET from towers}{with |ieta|=28}; MET (GeV)",met_bins,0,met_max);
    TL1RootHist* tower28MET=new TL1RootHist(tower28MET_template);
    plots.emplace(name,tower28MET);

    // Tower 28 MET Phi
    name = "tower28METPhi";
    TH1* tower28METPhi_template=new TH1F(name.c_str(),"#splitline{Phi of MET from towers}{with |ieta|=28}; Phi (#circ)",phi_bins,0,360);
    TL1RootHist* tower28METPhi=new TL1RootHist(tower28METPhi_template);
    plots.emplace(name,tower28METPhi);

    // MET in all towers
    name = "totalMET";
    TH1* totalMET_template=new TH1F(name.c_str(),"MET from all towers; MET (GeV)",met_bins,0,met_max);
    TL1RootHist* totalMET=new TL1RootHist(totalMET_template);
    plots.emplace(name,totalMET);

    // MET in all towers except 28
    name = "totalNot28MET";
    TH1* totalNot28MET_template=new TH1F(name.c_str(),"#splitline{MET from all towers}{except |ieta|=28}; MET (GeV)",met_bins,0,met_max);
    TL1RootHist* totalNot28MET=new TL1RootHist(totalNot28MET_template);
    plots.emplace(name,totalNot28MET);

    std::string overwrite_dir=dataset->baseOWdir+"/studyTower28MET/dists_";
    for(auto plot : plots){
        plot.second->SetOverwriteNames(overwrite_dir+outName+".root", name);
        plot.second->SetDrawOption("hist");
        plot.second->SetOutName(outName+"_"+plot.first);
    }
    return plots;
}

void FillPlots(int i_entry, PlotList& plots, const TL1EventClass* event){
        const int pu = event->GetPEvent()->fVertex->nVtx;
	//cout<<"Pile-up: "<<pu<<endl;
	// Skip events with less than 5 pile-up at this point
	if(pu<5) return; 

	const double& total_met=event->fRecalcL1Met;
	const double met_28=event->fMet28.met();
	const double met_28_phi=event->fMet28.phi()/TMath::Pi()*180;
	const double met_not_28=event->fMetNot28.met();

	auto iplot=plots.find("tower28MET");    if(iplot!=plots.end()) iplot->second->Fill(met_28,    1,pu);
         iplot=plots.find("totalMET");      if(iplot!=plots.end()) iplot->second->Fill(total_met, 1,pu);
         iplot=plots.find("totalNot28MET"); if(iplot!=plots.end()) iplot->second->Fill(met_not_28,1,pu);
    if (met_28>0){
         iplot=plots.find("tower28METPhi"); if(iplot!=plots.end()) iplot->second->Fill(met_28_phi,1,pu);
    }
}

#include "../MakePlots/makePlots.cxx"
MAIN_FUNCTION(studyTower28MET)
