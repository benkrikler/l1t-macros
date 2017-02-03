#include "../Plotting/TL1RootHists.h"
#include "../Plotting/TL1Rates.h"
#include "../Plotting/TL1Turnon.h"
#include "../Event/TL1EventClass.h"
#include "../Config/ntuple_cfg.h"
#include "../Config/sumTurnons_cfg.h"

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
    int met_bins=25;
    int phi_bins=36;

    // tmp variables that we use below
    std::string name;
    TH1* tmpl_hist;

    // Tower 28 MET
    name = "tower28METNotEmu";
    tmpl_hist=new TH1F(name.c_str(),"#splitline{MET from towers}{with |ieta|=28}; MET (GeV)",met_bins,0,met_max);
    plots.emplace(name,new TL1RootHist(tmpl_hist));

    // Tower 28 MET Emulated
    name = "tower28MET";
    tmpl_hist=new TH1F(name.c_str(),"#splitline{Emul. MET from towers}{with |ieta|=28}; MET (GeV)",met_bins,0,met_max);
    plots.emplace(name,new TL1RootHist(tmpl_hist));

    // Tower 28 MET Phi
    name = "tower28METPhi";
    tmpl_hist=new TH1F(name.c_str(),"#splitline{Phi of Emul. MET from towers}{with |ieta|=28}; Phi (#circ)",phi_bins,0,360);
    plots.emplace(name,new TL1RootHist(tmpl_hist));

    // MET in all towers
    name = "totalMETNotEmu";
    tmpl_hist=new TH1F(name.c_str(),"MET from all towers (not emulated); MET (GeV)",met_bins,0,met_max);
    plots.emplace(name,new TL1RootHist(tmpl_hist));

    // Emulated MET in all towers
    name = "totalMET";
    tmpl_hist=new TH1F(name.c_str(),"Emul. MET from all towers; MET (GeV)",met_bins,0,met_max);
    plots.emplace(name,new TL1RootHist(tmpl_hist));

    // MET in all towers except 28
    name = "totalNot28MET";
    tmpl_hist=new TH1F(name.c_str(),"#splitline{Emul. MET from all towers}{except |ieta|=28}; MET (GeV)",met_bins,0,met_max);
    plots.emplace(name,new TL1RootHist(tmpl_hist));

    std::string overwrite_dir=dataset->baseOWdir+"/studyTower28MET/dists_";
    for(auto plot : plots){
        plot.second->SetDrawOption("hist");
    }

    auto turnon_plots=sumTurnons(dataset);
    plots.insert(turnon_plots.begin(),turnon_plots.end());

    name = "turnonL1MetNo28";
    TL1Turnon* turnonL1MetNo28=new TL1Turnon;
    std::string seed = "recalcL1MetNot28";
    std::string xparam = "caloMetBERecalc";
    turnonL1MetNo28->SetSeed(seed, "L1 MET !28");
    turnonL1MetNo28->SetSeeds({0., 40., 60., 80., 100., 120.});
    turnonL1MetNo28->SetX(xparam, "Offline E_{T}^{miss} BE (GeV)");
    turnonL1MetNo28->SetXBins(metBins());
    turnonL1MetNo28->SetFit(dataset->doFit);
    plots.emplace(name,turnonL1MetNo28);

    for(auto plot : plots){
        plot.second->SetOverwriteNames(overwrite_dir+outName+".root", name);
        plot.second->SetOutName(outName+"_"+plot.first);
        plot.second->SetOutExtension("png");
    }

    return plots;
}

static int my_count=0, my_fill_count=0;

void FillPlots(int i_entry, PlotList& plots, const TL1EventClass* event, const ntuple_cfg* dataset){
    ++my_count;
    const int pu = event->GetPEvent()->fVertex->nVtx;
    if(pu<5) return; 

    if( dataset->triggerName == "SingleMu" and not event->fMuonFilterPassFlag ) return;
    if( not event->fMetFilterPassFlag ) return;
    ++my_fill_count;

    const double& total_met_notEmu=event->fRecalcL1Met;
    const double& total_met=event->fRecalcL1EmuMet;
    const double met_28_notEmu=event->fRecalcL1Met28.met();
    const double met_28=event->fRecalcL1MetEmu28.met();
    const double met_28_phi=event->fRecalcL1MetEmu28.phi()/TMath::Pi()*180;
    const double met_not_28=event->fRecalcL1MetEmuNot28.met();

    PlotList::iterator iplot;
    // MET distributions
    iplot=plots.find("tower28MET"     ); if(iplot!=plots.end()) iplot->second->Fill(met_28           , 1 , pu ); 
    iplot=plots.find("tower28METNotEmu");if(iplot!=plots.end()) iplot->second->Fill(met_28           , 1 , pu ); 
    iplot=plots.find("totalMETNotEmu" ); if(iplot!=plots.end()) iplot->second->Fill(total_met_notEmu , 1 , pu ); 
    iplot=plots.find("totalMET"       ); if(iplot!=plots.end()) iplot->second->Fill(total_met        , 1 , pu ); 
    iplot=plots.find("totalNot28MET"  ); if(iplot!=plots.end()) iplot->second->Fill(met_not_28       , 1 , pu ); 
    if (met_28>0){
        iplot=plots.find("tower28METPhi"); if(iplot!=plots.end()) iplot->second->Fill(met_28_phi     , 1 , pu );
    }

    // Turnon curves
    auto sums = event->GetPEvent()->fSums;
    const double& caloMetBE = sums->caloMetBE;
    iplot=plots.find("metBE"           ); if(iplot!=plots.end()) iplot->second->Fill(caloMetBE , event->fL1Met    , pu ); 
    iplot=plots.find("metBEEmu"        ); if(iplot!=plots.end()) iplot->second->Fill(caloMetBE , event->fL1EmuMet , pu ); 
    iplot=plots.find("metBERecalc"     ); if(iplot!=plots.end()) iplot->second->Fill(caloMetBE , total_met_notEmu , pu ); 
    iplot=plots.find("metBERecalcEmu"  ); if(iplot!=plots.end()) iplot->second->Fill(caloMetBE , total_met        , pu ); 
    iplot=plots.find("turnonL1MetNo28" ); if(iplot!=plots.end()) iplot->second->Fill(caloMetBE , met_not_28       , pu ); 
}

void Finalize(PlotList& plots){
    for(auto plot: plots){
        if(plot.first.find("Phi")!=0) TL1RootHist::SetLogY(false);
        else TL1RootHist::SetLogY(true);
        plot.second->DrawPlots();
        if(plot.first.find("met")==0 || plot.first.find("turnon")==0) continue;

        plot.second->NormaliseArea(1.);
        plot.second->SetDrawOption("histnostack");
        plot.second->DrawPlots("norm");
    }
    cout<<"Tried FillPlots() "<<my_count<<" times"<<endl;
    cout<<"Called Fill "<<my_fill_count<<" times"<<endl;
}

#include "../MakePlots/makePlots.cxx"
MAIN_FUNCTION(studyTower28MET)
