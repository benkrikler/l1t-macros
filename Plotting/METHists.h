#ifndef METHists_h
#define METHists_h


#include <string>

#include <TMarker.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TArrow.h>

#include "TL1Plots.h"
#include "TL1Turnon.h"
#include "TL1RootHists.h"
#include "../Debug/DebugHandler.h"
#include "../Config/ntuple_cfg.h"
#include "../Config/sumTurnons_cfg.h"

typedef std::map<std::string,TL1Plots*> PlotList;

class METHists:public TL1Plots{
    public:

        METHists(const std::string& name,const std::string& title,const std::string& leg):
                fName(name),fTitle(title),fLegendHeader(leg),fDataset(NULL),hTurnon(NULL),hDist(NULL){
            const double met_max=100;
            const int met_bins=25;

            TH1* tmpl_hist=new TH1F(("dist"+fName).c_str(),fTitle.c_str(),met_bins,0,met_max);
            hDist=new TL1RootHist(tmpl_hist);
            hDist->SetDrawOption("hist");

            hTurnon=new TL1Turnon;
            hTurnon->SetSeed(fName, fLegendHeader);
            hTurnon->SetSeeds({0.,70.,90., 110.});
            hTurnon->SetX("caloMetBERecalc", "Offline E_{T}^{miss} BE (GeV)");
            hTurnon->SetXBins(metBins());
        }
        void Register(const ntuple_cfg* dataset,PlotList& plots){
            fDataset=dataset;
            plots.emplace(fName,this);
        }

        virtual void InitPlots(){
            hTurnon->SetFit(fDataset->doFit);
            hTurnon->InitPlots();
            hDist->InitPlots();
        }
        virtual void OverwritePlots(){
            hTurnon->OverwritePlots();
            hDist->OverwritePlots();
        }
        virtual void Fill(const double & observedMet, const double & trueMet, const int & pu){
            if(!hTurnon) cout<<"hTurnon not initialised"<<endl;
            if(!hDist) cout<<"hDist not initialised"<<endl;
            if(!hTurnon || !hDist) return;

            hTurnon->Fill(trueMet,observedMet,pu);
            hDist->Fill(observedMet,1,pu);
        }
        virtual void DrawPlots(const char* name_append=NULL){
            hTurnon->DrawPlots(name_append);

            TL1RootHist::SetLogY(false);
            hDist->DrawPlots(name_append);
            hDist->NormaliseArea(1.);
            hDist->SetDrawOption("histnostack");
            std::string new_name_append="norm";
            if(name_append) new_name_append=name_append+new_name_append;
            hDist->DrawPlots(new_name_append.c_str());
        }
        virtual void SetOverwriteNames(const std::string & owRootName, const std::string & owHistName){
            if(hTurnon) hTurnon->SetOverwriteNames(owRootName,owHistName);
            if(hDist) hDist->SetOverwriteNames(owRootName,owHistName);}
        virtual void SetSample(const std::string & sampleName, const std::string & sampleTitle){
            if(hTurnon) hTurnon->SetSample(sampleName,sampleTitle);
            if(hDist) hDist->SetSample(sampleName,sampleTitle);}
        virtual void SetTrigger(const std::string & triggerName, const std::string & triggerTitle){
            if(hTurnon) hTurnon->SetTrigger(triggerName,triggerTitle);
            if(hDist) hDist->SetTrigger(triggerName,triggerTitle);}
        virtual void SetRun(const std::string & run){
            if(hTurnon) hTurnon->SetRun(run);
            if(hDist) hDist->SetRun(run);}
        virtual void SetOutName(const std::string & outName){
            if(hTurnon) hTurnon->SetOutName(outName);
            if(hDist) hDist->SetOutName(outName);}
        virtual void SetOutDir(const std::string & outDir){
            if(hTurnon) hTurnon->SetOutDir(outDir);
            if(hDist) hDist->SetOutDir(outDir);}
        virtual void SetOutExtension(const std::string & outExt){
            if(hTurnon) hTurnon->SetOutExtension(outExt);
            if(hDist) hDist->SetOutExtension(outExt);}
        virtual void SetAddMark(const std::string & addMark){
            if(hTurnon) hTurnon->SetAddMark(addMark);
            if(hDist) hDist->SetAddMark(addMark);}
        virtual void SetPuType(const std::vector<std::string> & puType){
            if(hTurnon) hTurnon->SetPuType(puType);
            if(hDist) hDist->SetPuType(puType);}
        virtual void SetPuBins(const std::vector<int> & puBins){
            if(hTurnon) hTurnon->SetPuBins(puBins);
            if(hDist) hDist->SetPuBins(puBins);}
        virtual void SetPuFile(const std::string & puFileName){
            if(hTurnon) hTurnon->SetPuFile(puFileName);
            if(hDist) hDist->SetPuFile(puFileName);}

    public:
        std::string fName, fTitle,fLegendHeader;
        const ntuple_cfg* fDataset;
        TL1Turnon* hTurnon;
        TL1RootHist* hDist;
};

#endif // METHists_h
