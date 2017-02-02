#ifndef TL1RootHist_h
#define TL1RootHist_h

#include <string>
#include <vector>
#include <sstream>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TLatex.h>
#include <TStyle.h>

#include "TL1Plots.h"
#include "../Debug/DebugHandler.h"

class TL1RootHist : public TL1Plots
{
    public:
        TL1RootHist():
		TL1Plots(),fTemplatePlot(NULL),fRootFile(NULL),fNDimensions(0){}
        TL1RootHist(TH1* tmpl):
		TL1Plots(),fTemplatePlot(NULL),fRootFile(NULL),fNDimensions(0){
		SetTemplatePlot(tmpl);
	}
        ~TL1RootHist();
        
        virtual void InitPlots();
        virtual void OverwritePlots();
        virtual void Fill(const double & xVal, const double & yVal, const int & pu);
        virtual void DrawPlots();
        void DrawCmsStamp();
        void SetTemplatePlot(TH1* tmplt);
        TH1* GetTemplatePlot()const{return fTemplatePlot;}
        int GetNDimensions()const{return fNDimensions;}

    private:
        void MakePlot(const TString name);
        void FillPlot(TH1* hist,const double xVal,const double yVal,const double pu_weight);
    private:
        TH1* fTemplatePlot;
        std::vector<TH1*> fPlot;
        TFile * fRootFile;
        int fNDimensions;

};

TL1RootHist::~TL1RootHist()
{
    fRootFile->Close();
    delete fRootFile;
}

void TL1RootHist::SetTemplatePlot(TH1* tmplt){
    fTemplatePlot=tmplt;
    fNDimensions=0;
         if(fTemplatePlot->InheritsFrom(TH3::Class())) fNDimensions=3;
    else if(fTemplatePlot->InheritsFrom(TH2::Class())) fNDimensions=2;
    else if(fTemplatePlot->InheritsFrom(TH1::Class())) fNDimensions=1;
}

void TL1RootHist::MakePlot(const TString name){
    fPlot.push_back((TH1*)fTemplatePlot->Clone(Form("%s_%s",fTemplatePlot->GetName(),name.Data())));
    fPlot.back()->Sumw2();
    fPlot.back()->SetDirectory(0);
}

void TL1RootHist::InitPlots() {
    fRootFile = TFile::Open(Form("%s/%s.root", this->GetOutDir().c_str(), this->GetOutName().c_str()), "RECREATE");

    MakePlot("all");
    for(unsigned int ipu=0; ipu<this->GetPuType().size(); ++ipu) {
        MakePlot(this->GetPuType()[ipu]);
    }
}

void TL1RootHist::OverwritePlots()
{
    fPlot.clear();
    TFile * rootFile = TFile::Open(this->GetOverwriteRootFilename().c_str(), "READ");
    fRootFile = TFile::Open(Form("%s/%s_overwrite.root",this->GetOutDir().c_str(),this->GetOutName().c_str()),"RECREATE");

    fPlot.push_back((TH1F*)rootFile->Get(this->GetOverwriteHistname().c_str()));
    for(unsigned int ipu=0; ipu<this->GetPuType().size(); ++ipu)
    {
        fPlot.push_back((TH1F*)rootFile->Get(Form("%s_%s",this->GetOverwriteHistname().c_str(),this->GetPuType()[ipu].c_str())));
        fPlot.back()->SetDirectory(0);
    }
    rootFile->Close();
    delete rootFile;
}

void TL1RootHist::FillPlot(TH1* hist,const double xVal,const double yVal,const double pu_weight){
   switch (fNDimensions){
       case 1:
          hist->Fill(xVal,pu_weight);
            break;
      case 2: case 3:
          static_cast<TH2*>(hist)->Fill(xVal,yVal,pu_weight);
          break;
      default:
          cout<<"TL1RootHist::FillPlot() Error: bad number of dimensions to fill plot, "<<hist->GetName()<<endl;
  }
}

void TL1RootHist::Fill(const double & xVal, const double & yVal, const int & pu) {
    const double pu_weight=this->GetPuWeight(pu);
    FillPlot(fPlot[0],xVal,yVal,pu_weight);
    for(unsigned int ipu=0; ipu<this->GetPuType().size(); ++ipu) {
        if( pu >= this->GetPuBins()[ipu] && pu < this->GetPuBins()[ipu+1] ){
            FillPlot(fPlot[ipu],xVal,yVal,pu_weight);
            break;
        }
    }
}

void TL1RootHist::DrawPlots()
{
    TCanvas * can(new TCanvas(Form("can_%f",this->GetRnd()),""));
    
    fRootFile->WriteTObject(fPlot[0]);
    fPlot[0]->Draw();
    can->SetLogy();
    DrawCmsStamp();
    TLatex* title=new TLatex(0.64,0.87,fTemplatePlot->GetTitle());
    title->SetNDC();
    title->SetTextAlign(33);
    title->Draw();

    std::string outName = Form("%s/%s.%s",this->GetOutDir().c_str(),this->GetOutName().c_str(),this->GetOutExtension().c_str());
    can->SaveAs(outName.c_str());

    if( this->GetPuType().size() <= 0 ) return;

    TCanvas * can2(new TCanvas(Form("can_%f",this->GetRnd()),""));
    TLegend * leg2(new TLegend(0.65,0.55,0.88,0.55+0.05*this->GetPuType().size()));
    THStack* stack=new THStack;
    for(unsigned int ipu=0; ipu<this->GetPuType().size(); ++ipu)
    {
        TH1* plot=fPlot[ipu+1];
        fRootFile->WriteTObject(plot);
        std::stringstream entryName;
        entryName << this->GetPuBins()[ipu] << " #leq PU";
        if( ipu<this->GetPuType().size()-1 ) entryName<<" < " << this->GetPuBins()[ipu+1];
        SetColor(plot,ipu, this->GetPuType().size(),true);
        leg2->AddEntry(plot, entryName.str().c_str());
        stack->Add(plot);
    }
    stack->Draw(GetDrawOption().c_str());
    stack->GetXaxis()->SetTitle(fTemplatePlot->GetXaxis()->GetTitle());
    stack->GetYaxis()->SetTitle(fTemplatePlot->GetYaxis()->GetTitle());
    can2->SetLogy();
    DrawCmsStamp();
    leg2->Draw();
    title->DrawClone();
    can2->Update();

    outName = Form("%s/%s_puBins.%s", this->GetOutDir().c_str(), this->GetOutName().c_str(),this->GetOutExtension().c_str());
    can2->SaveAs(outName.c_str());
    delete can;
    delete can2;
}

void TL1RootHist::DrawCmsStamp()
{
    TLatex * latex(new TLatex());
    latex->SetNDC();
    latex->SetTextFont(42);
    if( this->GetSampleName() == "Data" )
    {
        latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Preliminary} 2016 Data");
        latex->SetTextAlign(31);
        latex->DrawLatex(0.92,0.92,Form("%s (13 TeV)",this->GetRun().c_str()));
    }
    else
    {
        latex->DrawLatex(0.15,0.92,"#bf{CMS} #it{Simulation Preliminary}");
        latex->SetTextAlign(31);
        latex->DrawLatex(0.92,0.92,Form("%s (13 TeV)",this->GetSampleName().c_str()));
    }
    latex->SetTextAlign(32);
    latex->DrawLatex(0.87,0.82,this->GetAddMark().c_str());
}

#endif // TL1RootHist_h
