#include "../Plotting/TL1RootHists.h"
#include "../Plotting/TL1Rates.h"
#include "../Plotting/TL1Turnon.h"
#include "../Plotting/METHists.h"
#include "../Event/TL1EventClass.h"
#include "../Config/ntuple_cfg.h"
#include "../Config/sumTurnons_cfg.h"

#include <string>
#include <iostream>
#include <TH1F.h>

using std::cout;
using std::endl;

//void MakePlots(ntuple_cfg*,PlotList& plots);
void FillPlots(int i_entry, PlotList& plots, const TL1EventClass*);

void MakePlots(ntuple_cfg* dataset,PlotList& plots){
    plots.clear();

    METHists* recalcL1EmuMet           =new METHists("RecalcL1EmuMet"           , "Default MET, all barrel towers"                   , "L1MET(Barrel)"        ); 
    METHists* recalcL1EmuMetHF         =new METHists("RecalcL1EmuMetHF"         , "Default MET, all barrel + forward towers"         , "L1MET(Barrel+Fwd)"    ); 
    METHists* recalcL1EmuMet28Only     =new METHists("RecalcL1EmuMet28Only"     , "MET in only |ieta|=28"                            , "L1MET(|ieta|=28)"    ); 
    METHists* recalcL1EmuMetNot28      =new METHists("RecalcL1EmuMetNot28"      , "MET with |ieta|<28"                               , "L1MET(|ieta|<28)"  ); 
    METHists* recalcL1EmuMetPUS        =new METHists("RecalcL1EmuMetPUS"        , "MET with PUS"                                     , "L1MET(Barrel) PUS"    ); 
    METHists* recalcL1EmuMetPUSHF      =new METHists("RecalcL1EmuMetPUSHF"      , "MET with PUS + forward towers"                    , "L1MET(Barrel+Fwd) PUS"); 
    METHists* recalcL1EmuMetPUS28      =new METHists("RecalcL1EmuMetPUS28"      , "MET with PUS for |ieta|=28"                       , "L1MET(Barrel) PUS28"  ); 
    METHists* recalcL1EmuMetPUSThresh  =new METHists("RecalcL1EmuMetPUSThresh"  , "MET with PUS & E_{T} threshold"                   , "L1MET(Barrel) PUS+Threshold"); 
    METHists* recalcL1EmuMetPUSThreshHF=new METHists("RecalcL1EmuMetPUSThreshHF", "MET with PUS & E_{T} threshold and forward towers", "L1MET(Barrel+Fwd) PUS28+Threshold"  ); 
    METHists* recalcL1Met              =new METHists("RecalcL1Met"           , "(Not Emu.) Default MET, all barrel towers"              , "L1MET raw (Barrel)"        ); 
    METHists* recalcL1Met28Only        =new METHists("RecalcL1Met28Only"     , "(Not Emu.) MET in only |ieta|=28"                       , "L1MET raw (|ieta|=28)"    ); 

    recalcL1EmuMet           ->Register(dataset,plots);
    recalcL1EmuMetHF         ->Register(dataset,plots);
    recalcL1EmuMet28Only     ->Register(dataset,plots);
    recalcL1EmuMetNot28      ->Register(dataset,plots);
    recalcL1EmuMetPUS        ->Register(dataset,plots);
    recalcL1EmuMetPUSHF      ->Register(dataset,plots);
    recalcL1EmuMetPUS28      ->Register(dataset,plots);
    recalcL1EmuMetPUSThresh  ->Register(dataset,plots);
    recalcL1EmuMetPUSThreshHF->Register(dataset,plots); 
    recalcL1Met              ->Register(dataset,plots);
    recalcL1Met28Only        ->Register(dataset,plots);

    // Tower 28 MET Phi
    int phi_bins=18;
    std::string name = "tower28METPhi";
    TH1* tmpl_hist=new TH1F(name.c_str(),"#splitline{Phi of Emul. MET from towers}{with |ieta|=28}; Phi (#circ)",phi_bins,0,360);
    plots.emplace(name,new TL1RootHist(tmpl_hist));

    auto turnon_plots=sumTurnons(dataset);
    plots.insert(turnon_plots.begin(),turnon_plots.end());

    std::string overwrite_dir=dataset->baseOWdir+"/studyTower28MET/";
    std::string outName = dataset->triggerName;
    for(auto plot : plots){
        plot.second->SetOverwriteNames(overwrite_dir+outName+".root", name);
        plot.second->SetOutName(outName+"_"+plot.first);
        plot.second->SetOutExtension("pdf");
    }
}

static int my_count=0, my_fill_count=0, fRecalcL1Met_nonZero=0;

void FillPlots(int i_entry, PlotList& plots, const TL1EventClass* event, const ntuple_cfg* dataset){
    ++my_count;
    const int pu = event->GetPEvent()->fVertex->nVtx;
    if(pu<5) return; 

    if( dataset->triggerName == "SingleMu" and not event->fMuonFilterPassFlag ) return;
    if( not event->fMetFilterPassFlag ) return;
    ++my_fill_count;
    if(event->fRecalcL1Met >0) ++fRecalcL1Met_nonZero;

    auto sums = event->GetPEvent()->fSums;
    const double& caloMetBE = sums->caloMetBE;
    PlotList::iterator iplot;
    iplot=plots.find("RecalcL1EmuMet"           ); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1EmuMet.met()           ,caloMetBE, pu ); 
    iplot=plots.find("RecalcL1EmuMetHF"         ); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1EmuMetHF.met()         ,caloMetBE, pu ); 
    iplot=plots.find("RecalcL1EmuMet28Only"     ); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1EmuMet28Only.met()     ,caloMetBE, pu ); 
    iplot=plots.find("RecalcL1EmuMetNot28"      ); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1EmuMetNot28.met()      ,caloMetBE, pu ); 
    iplot=plots.find("RecalcL1EmuMetPUS"        ); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1EmuMetPUS.met()        ,caloMetBE, pu ); 
    iplot=plots.find("RecalcL1EmuMetPUSHF"      ); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1EmuMetPUSHF.met()      ,caloMetBE, pu ); 
    iplot=plots.find("RecalcL1EmuMetPUS28"      ); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1EmuMetPUS28.met()      ,caloMetBE, pu ); 
    iplot=plots.find("RecalcL1EmuMetPUSThresh"  ); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1EmuMetPUSThresh.met()  ,caloMetBE, pu ); 
    iplot=plots.find("RecalcL1EmuMetPUSThreshHF"); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1EmuMetPUSThreshHF.met(),caloMetBE, pu ); 
    iplot=plots.find("RecalcL1Met"              ); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1Met                    ,caloMetBE, pu ); 
    iplot=plots.find("RecalcL1Met28Only"        ); if(iplot!=plots.end()) iplot->second->Fill(event->fRecalcL1Met28Only.met()        ,caloMetBE, pu ); 

    const double met_28=event->fRecalcL1EmuMet28Only.met();
    if (met_28>0){
        const double met_28_phi=event->fRecalcL1EmuMet28Only.phi()/TMath::Pi()*180;
        iplot=plots.find("tower28METPhi"); if(iplot!=plots.end()) iplot->second->Fill(met_28_phi     , 1 , pu );
    }

    // Turnon curves
    iplot=plots.find("metBE"           ); if(iplot!=plots.end()) iplot->second->Fill(caloMetBE , event->fL1Met                , pu );
    iplot=plots.find("metBEEmu"        ); if(iplot!=plots.end()) iplot->second->Fill(caloMetBE , event->fL1EmuMet             , pu );
    iplot=plots.find("metBERecalc"     ); if(iplot!=plots.end()) iplot->second->Fill(caloMetBE , event->fRecalcL1Met          , pu );
    iplot=plots.find("metBERecalcEmu"  ); if(iplot!=plots.end()) iplot->second->Fill(caloMetBE , event->fRecalcL1EmuMet.met() , pu );
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
    cout<<"fRecalcL1Met_nonZero:  "<<fRecalcL1Met_nonZero<<endl;
}

#include "../MakePlots/makePlots.cxx"
MAIN_FUNCTION(studyTower28MET)
