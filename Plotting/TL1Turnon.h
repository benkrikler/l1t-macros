#ifndef TL1TURNON_H
#define TL1TURNON_H

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
#include "../Debug/DebugHandler.h"

TGraphAsymmErrors GetEfficiency(TH1F * total, TH1F * pass);

class TL1Turnon : public TL1Plots
{
    public:
        ~TL1Turnon();

        virtual void InitPlots();
        virtual void OverwritePlots();
        virtual void Fill(const double & xVal, const double & seedVal, const int & pu);
        virtual void DrawPlots(const char* name_append=NULL);

        void SetSeeds(const vector<double> & seeds){ fSeeds = seeds; }
        void SetXBins(const vector<double> & xBins){ fXBins = xBins; }
        void SetX(const std::string & xName, const std::string & xTitle){ fXName = xName; fXTitle = xTitle; }
        void SetSeed(const std::string & seedName, const std::string & seedTitle){ fSeedName = seedName; fSeedTitle = seedTitle; }
        void SetFit(const bool & doFit){ fDoFit = doFit; }

    private:
        void DrawCmsStamp(std::string stampPos="Left");
        void DrawTurnons();
        void DrawFitResults();
        void DrawCmsStampTurnon(const double & max);
        TF1* fit(TGraphAsymmErrors * eff, double p50);

    private:
        std::vector<std::vector<TH1F*>> fPlots;
        std::vector<std::vector<TGraphAsymmErrors*>> fTurnons;
        std::vector<std::vector<TF1*>> fFits;

        TFile* fPlotsRoot;
        TFile* fTurnonsRoot;

        vector<double> fSeeds, fXBins;
        std::string fXName, fSeedName;
        std::string fXTitle, fSeedTitle;
        bool fDoFit;
};

TL1Turnon::~TL1Turnon() {
    delete fPlotsRoot;
    delete fTurnonsRoot;
}

void TL1Turnon::InitPlots() {
    fPlotsRoot = TFile::Open(Form("%s/dists_%s.root",this->GetOutDir().c_str(),this->GetOutName().c_str()),"RECREATE");
    fTurnonsRoot = TFile::Open(Form("%s/effs_%s.root",this->GetOutDir().c_str(),this->GetOutName().c_str()),"RECREATE");

    for(unsigned i=0; i<fSeeds.size(); ++i)
    {
        std::vector<TH1F*> temp;
        temp.emplace_back(new TH1F(Form("dist_%s_%s_%g",fXName.c_str(),fSeedName.c_str(),fSeeds[i]),"", fXBins.size()-1,&(fXBins)[0]));
        temp.back()->SetDirectory(0);
        temp.back()->Sumw2();
        temp.back()->GetXaxis()->SetTitle(fXTitle.c_str());
        temp.back()->GetYaxis()->SetTitle("Number of Entries");
        this->SetColor(temp.back(), i-1, fSeeds.size()-1);

        for(int ipu=0; ipu<this->GetPuType().size(); ++ipu)
        {
            temp.emplace_back(new TH1F(Form("dist_%s_%s_%g_%s",fXName.c_str(),fSeedName.c_str(),fSeeds[i],this->GetPuType()[ipu].c_str()),"", fXBins.size()-1,&(fXBins)[0]));
            temp.back()->SetDirectory(0);
            temp.back()->Sumw2();
            temp.back()->GetXaxis()->SetTitle(fXTitle.c_str());
            temp.back()->GetYaxis()->SetTitle("Number of Entries");
            this->SetColor(temp.back(), ipu, this->GetPuType().size());
        }
        fPlots.push_back(temp);
    }
}

void TL1Turnon::OverwritePlots()
{
    fPlots.clear();
    TFile * rootFile = TFile::Open(this->GetOverwriteRootFilename().c_str(),"READ");

    fPlotsRoot = TFile::Open(Form("%s/dists_%s_overwrite.root",this->GetOutDir().c_str(),this->GetOutName().c_str()),"RECREATE");
    fTurnonsRoot = TFile::Open(Form("%s/effs_%s_overwrite.root",this->GetOutDir().c_str(),this->GetOutName().c_str()),"RECREATE");

    for(unsigned i=0; i<fSeeds.size(); ++i)
    {
        //cout << "fSeeds[i] = " << fSeeds[i] << endl;
        std::vector<TH1F*> temp;
        temp.push_back((TH1F*)rootFile->Get(Form("%s_%i",this->GetOverwriteHistname().c_str(),(int)fSeeds[i])));
        temp.back()->SetDirectory(0);
        temp.back()->GetXaxis()->SetTitle(fXTitle.c_str());
        temp.back()->GetYaxis()->SetTitle("Number of Entries");
        this->SetColor(temp.back(), i-1, fSeeds.size()-1);

        for(int ipu=0; ipu<this->GetPuType().size(); ++ipu)
        {
            //cout << "GetPuType()[ipu] = " << this->GetPuType()[ipu] << endl;
            temp.push_back((TH1F*)rootFile->Get(Form("%s_%i_%s",this->GetOverwriteHistname().c_str(),(int)fSeeds[i],this->GetPuType()[ipu].c_str())));
            temp.back()->SetDirectory(0);
            temp.back()->GetXaxis()->SetTitle(fXTitle.c_str());
            temp.back()->GetYaxis()->SetTitle("Number of Entries");
            this->SetColor(temp.back(), ipu, this->GetPuType().size());
        }
        fPlots.push_back(temp);
    }
    delete rootFile;
}

void TL1Turnon::Fill(const double & xVal, const double & seedVal, const int & pu)
{
    for(unsigned i=0; i<fSeeds.size(); ++i)
    {
        if( !(seedVal >= fSeeds[i]) ) break;
        fPlots[i][0]->Fill(xVal,this->GetPuWeight(pu));

        for(int ipu=0; ipu<this->GetPuType().size(); ++ipu)
        {
            if( pu >= this->GetPuBins()[ipu] && pu < this->GetPuBins()[ipu+1] )
                fPlots[i][ipu+1]->Fill(xVal,this->GetPuWeight(pu));
        }
    }
}

void TL1Turnon::DrawPlots(const char* name_append)
{
    TCanvas * can(new TCanvas(Form("can_%f",this->GetRnd()),"")); 
    TLegend * leg(new TLegend(0.58,0.35,0.88,0.55));
    for(unsigned i=0; i<fPlots.size(); ++i)
    {
        if(i==0) fPlots[i][0]->Draw();
        else fPlots[i][0]->Draw("same");
        fPlotsRoot->WriteTObject(fPlots[i][0]);
        leg->AddEntry(fPlots[i][0], Form("%s > %g GeV", fSeedTitle.c_str(), fSeeds[i]));

        for(int ipu=0; ipu<this->GetPuType().size(); ++ipu)
        {
            fPlots[i][ipu+1]->Draw("same");
            fPlotsRoot->WriteTObject(fPlots[i][ipu+1]);

            std::stringstream entryName;
            if( ipu < this->GetPuType().size()-1 ) entryName << this->GetPuBins()[ipu] << " #leq PU < " << this->GetPuBins()[ipu+1];
            else entryName << this->GetPuBins()[ipu] << " #leq PU";
            leg->AddEntry(fPlots[i][ipu+1],entryName.str().c_str());
            entryName.str("");
        }
    }
    can->SetLogy();
    leg->Draw();

    DrawCmsStamp();

    std::string outName = Form("%s/dists_%s.%s",this->GetOutDir().c_str(),this->GetOutName().c_str(),this->GetOutExtension().c_str());
    can->SaveAs(outName.c_str());
    delete can;

    cout<<"Drawing turnons: "<<endl;
    DrawTurnons();
    if(fDoFit){
        cout<<"Drawing fit results: "<<endl;
        DrawFitResults();
    }
}

void TL1Turnon::DrawCmsStamp(std::string stampPos)
{
    TLatex * latex(new TLatex());
    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextAlign(32);
    latex->DrawLatexNDC(0.18,0.92,"#bf{CMS} #it{Preliminary}");
    latex->DrawLatexNDC(0.92,0.92,"(13 TeV)");
    if( this->GetSampleName() == "Data" )
    {
        //latex->DrawLatex(0.89,0.80,"#it{Preliminary}");
        latex->SetTextAlign(31);
        std::string runNo = "run " + this->GetRun() + ", ";
        //latex->DrawLatex(0.92, 0.92, Form("%s%s, #sqrt{s} = 13 TeV",runNo.c_str(),this->GetTriggerTitle().c_str()));
        //latex->DrawLatex(0.92,0.92,"(13 TeV)");
    }
    else
    {
        latex->DrawLatexNDC(0.89,0.80,"#it{Simulation}");
        latex->DrawLatexNDC(0.89,0.75,"#it{Preliminary}");
        latex->SetTextAlign(31);
        latex->DrawLatexNDC(0.92, 0.92, Form("%s, #sqrt{s} = 13 TeV",this->GetSampleTitle().c_str()));
    }
    latex->SetTextAlign(11);
    //latex->DrawLatex(0.18,0.92,this->GetAddMark().c_str());
}

void TL1Turnon::DrawTurnons() {
    fFits.clear();
    TCanvas * nomCan(new TCanvas(Form("can_%f",this->GetRnd()),"c1"));
    TLegend * nomLeg(new TLegend(0.58,0.15,0.83,0.15+0.16*(2+fSeeds.size())/5.0,fSeedTitle.c_str()));
    nomLeg->SetTextSize(0.04);
    TArrow * arrow = new TArrow();
    double max(0.0);
    for(int i=1; i<fSeeds.size(); ++i)
    {
        std::vector<TGraphAsymmErrors*> temp;
        temp.emplace_back(new TGraphAsymmErrors(GetEfficiency(fPlots[0][0], fPlots[i][0])));
        arrow->SetLineColor(fPlots[i][0]->GetLineColor());
        arrow->SetFillColor(fPlots[i][0]->GetLineColor());
        temp[0]->SetLineColor(fPlots[i][0]->GetLineColor());
        temp[0]->SetMarkerColor(fPlots[i][0]->GetMarkerColor());
        temp[0]->SetFillColor(0);
        temp[0]->GetXaxis()->SetTitle(fPlots[i][0]->GetXaxis()->GetTitle());
        max = temp[0]->GetX()[temp[0]->GetN()-1]+0.9*temp[0]->GetErrorXhigh(temp[0]->GetN()-1);
        temp[0]->GetXaxis()->SetLimits(temp[0]->GetX()[0]-temp[0]->GetErrorXlow(0), max);
        temp[0]->GetYaxis()->SetTitle("Efficiency");
        temp[0]->SetMinimum(0.0);
        temp[0]->SetMaximum(1.1);
        nomCan->cd();
        //temp[0]->GetXaxis()->SetRangeUser(100,3500);
        arrow->DrawArrow(temp[0]->GetX()[temp[0]->GetN()-1]+0.89*temp[0]->GetErrorXhigh(temp[0]->GetN()-1),
                temp[0]->GetY()[temp[0]->GetN()-1], max,
                temp[0]->GetY()[temp[0]->GetN()-1],
                0.013);
        fTurnonsRoot->WriteTObject(temp[0]);

        if( fDoFit ){
            TF1* tmp_fit=fit(temp[0], fSeeds[i]);
            //tmp_fit->Draw("apsame");
            fFits.emplace_back();
            fFits.back().emplace_back(tmp_fit);
            fTurnonsRoot->WriteTObject(tmp_fit);
        }
        if( i == 1 ) temp[0]->Draw("ap");
        else temp[0]->Draw("psame");
        nomLeg->AddEntry(temp[0], Form("%g",fSeeds[i]));

        TCanvas * puCan(new TCanvas(Form("puCan_%f",this->GetRnd()),""));
        TLegend * puLeg(new TLegend(0.65,0.15,0.9,0.15+0.08*this->GetPuType().size(),Form("%s > %g",fSeedTitle.c_str(),fSeeds[i])));
        puLeg->SetTextSize(0.04);
        double puMax(0.0);
        for(int ipu=0; ipu<GetPuType().size(); ++ipu)
        {
            temp.emplace_back(new TGraphAsymmErrors(GetEfficiency(fPlots[0][ipu+1], fPlots[i][ipu+1])));
            arrow->SetLineColor(fPlots[i][ipu+1]->GetLineColor());
            arrow->SetFillColor(fPlots[i][ipu+1]->GetLineColor());
            temp[ipu+1]->SetLineColor(fPlots[i][ipu+1]->GetLineColor());
            temp[ipu+1]->SetMarkerColor(fPlots[i][ipu+1]->GetMarkerColor());
            temp[ipu+1]->SetFillColor(0);
            temp[ipu+1]->GetXaxis()->SetTitle(fPlots[i][ipu+1]->GetXaxis()->GetTitle());
            puMax = temp[ipu+1]->GetX()[temp[ipu+1]->GetN()-1]+0.9*temp[ipu+1]->GetErrorXhigh(temp[ipu+1]->GetN()-1);
            temp[ipu+1]->GetXaxis()->SetLimits(temp[ipu+1]->GetX()[0]-temp[ipu+1]->GetErrorXlow(0), puMax);
            temp[ipu+1]->GetYaxis()->SetTitle("Efficiency");
            temp[ipu+1]->SetMinimum(0.0);
            temp[ipu+1]->SetMaximum(1.1);
            puCan->cd();
            arrow->DrawArrow(temp[ipu]->GetX()[temp[ipu]->GetN()-1]+0.89*temp[ipu]->GetErrorXhigh(temp[ipu]->GetN()-1),
                    temp[ipu]->GetY()[temp[ipu]->GetN()-1], puMax,
                    temp[ipu]->GetY()[temp[ipu]->GetN()-1],
                    0.013);
            fTurnonsRoot->WriteTObject(temp[ipu+1]);

            if( fDoFit ) {
                TF1* tmp_fit=fit(temp[ipu+1], fSeeds[i]);
                //tmp_fit->Draw("apsame");
                fTurnonsRoot->WriteTObject(tmp_fit);
                fFits.back().emplace_back(tmp_fit);
            }
            if( ipu == 0 ) temp[ipu+1]->Draw("ap");
            else temp[ipu+1]->Draw("psame");

            std::stringstream entryName;
            if( ipu < this->GetPuType().size()-1 ) entryName << this->GetPuBins()[ipu] << " #leq PU < " << this->GetPuBins()[ipu+1];
            else entryName << this->GetPuBins()[ipu] << " #leq PU";
            puLeg->AddEntry(temp[ipu+1],entryName.str().c_str());
            entryName.str("");
        }
        puLeg->Draw();
        DrawCmsStampTurnon(puMax);
        TLatex * latex = new TLatex();
        latex->SetNDC();
        //latex->SetTextFont(42);
        //latex->DrawLatex(0.65,0.41,this->GetAddMark().c_str());

        std::string puOutName = Form("%s/effs_%s_puBins_seed%i.%s",this->GetOutDir().c_str(),this->GetOutName().c_str(),(int)fSeeds[i],this->GetOutExtension().c_str());
        puCan->SaveAs(puOutName.c_str());
        delete puCan;
    }
    nomCan->cd();
    nomLeg->Draw();
    DrawCmsStampTurnon(max);

    TLatex * nomlatex = new TLatex();
    nomlatex->SetNDC();
    nomlatex->SetTextFont(42);
    nomlatex->SetTextAlign(31);
    //nomlatex->DrawLatex(0.8,0.15+0.2*(2+fSeeds.size())/5.0+0.02,"<PU>=14");

    std::string nomOutName = Form("%s/effs_%s.%s",this->GetOutDir().c_str(),this->GetOutName().c_str(),this->GetOutExtension().c_str());
    nomCan->SaveAs(nomOutName.c_str());
    delete nomCan;
}

void TL1Turnon::DrawFitResults(){
    TCanvas * can(new TCanvas(Form("can_%f",this->GetRnd()),"c1"));
    TGraphErrors* errors[fSeeds.size()-1];
    TGraphErrors* errors_allPu[fSeeds.size()-1];
    TLegend * leg  =   new TLegend(0.17,0.72,0.9,0.78);
    TLegend * seed_leg=new TLegend(0.17,0.78,0.9,0.84);
    leg->SetNColumns(3);
    seed_leg->SetNColumns(fSeeds.size()-1);
    seed_leg->SetMargin(0.5);
    leg->SetFillStyle(1001);
    seed_leg->SetFillStyle(1001);
    leg->SetTextSize(0.04);
    seed_leg->SetTextSize(0.04);

    for(int i_seed=0; i_seed<fSeeds.size()-1; ++i_seed){
        errors[i_seed]=new TGraphErrors(GetPuType().size());
        errors_allPu[i_seed]=new TGraphErrors(1);
        errors_allPu[i_seed]->SetMarkerStyle(kOpenCircle);
        //errors_allPu[i_seed]->SetMarkerSize(2);
        SetColor(errors[i_seed],i_seed,fSeeds.size()-0.5,true);
        //SetColor(errors_allPu[i_seed],i_seed+0.5,fSeeds.size()-0.5,false);
        seed_leg->AddEntry(errors[i_seed],Form("%g",fSeeds[i_seed+1]),"lep");

        TMarker* marker_start=new TMarker(0,0,kOpenSquare);
        TMarker* marker_stop =new TMarker(0,0,kOpenTriangleUp);
        marker_start->SetMarkerSize(2);
        marker_stop ->SetMarkerSize(2);
        errors[i_seed]->GetListOfFunctions()->Add(marker_start);
        errors[i_seed]->GetListOfFunctions()->Add(marker_stop );
        if(i_seed==0){
            //leg->AddEntry((TObject*)NULL,fSeedTitle.c_str());
            leg->AddEntry(marker_start,"First Pile-up Bin","p");
            leg->AddEntry(errors_allPu[i_seed],"All Pile-up","p");
            leg->AddEntry(marker_stop,"Last Pile-up Bin","p");
            TText* title=new TLatex(0.17,0.88,fSeedTitle.c_str());
            title->SetNDC();
            title->SetTextAlign(13);
            title->SetTextFont(72);
            title->SetTextSize(0.05);
            errors[i_seed]->GetListOfFunctions()->Add(title);
        }

        int point=0;
        for(int ipu=0; ipu<GetPuType().size()+1; ++ipu){
            TF1* params=fFits[i_seed][ipu];
            const int n_pars=4;
            double par[n_pars],err[n_pars];
            for (int i=0;i<n_pars; ++i) {
                par[i]=params->GetParameter(i);
                err[i]=params->GetParError(i);
            }
            const double mu=par[0];
            const double asymm=par[2];
            const double mu_err=err[0];
            const double asymm_err=err[2];
            //const double sigma=1/par[1];
            //const double sigma_err=err[1]/sigma/sigma;
            const double sigma=TMath::Power(par[1]*asymm,-0.5);
            const double sigma_err=0.5*sigma*sigma*sigma*(asymm*err[1] + par[1]*asymm_err);
            //const double sigma=TMath::Sqrt(par[1]/asymm);
            //const double sigma_err=0.5/asymm/sigma*(err[1] + sigma*sigma*asymm_err);

            //const double sigma=1/par[0];
            //const double mu=par[1];
            //const double asymm=par[2];
            //const double sigma_err=sigma*sigma*err[0];
            //const double mu_err=err[1];
            //const double asymm_err=err[2];
            cout<<"fittable: "<<GetOutName()<<": "<<i_seed<<" ("<<fSeeds[i_seed+1]<<") "<<ipu<<" "<<mu<<"+-"<<mu_err<<"  "<<sigma<<"+-"<<sigma_err<<"  "<<asymm<<"+-"<<asymm_err<<" "<<par[3]<<"+-"<<err[3]<<endl;
            TGraphErrors* plot=NULL;
            if(ipu==0) plot=errors_allPu[i_seed];
            else {
                plot=errors[i_seed];
                if(mu!=mu or sigma!=sigma ) continue; // NaN checks
                if(mu==0 and sigma==0 ) continue;
                if(sigma > 60 ) continue;
                if(mu_err > 5 ) continue;
                if( fabs(fSeeds[i_seed+1]-mu)>30) continue;
            }

            plot->SetPoint(point,mu,sigma);
            plot->SetPointError(point,mu_err,sigma_err);
            if(ipu==0) continue; // skip the next steps if this is not in a puBin
            if(point==0) {
                marker_start->SetX(mu); marker_start->SetY(sigma); 
            }
            marker_stop->SetX(mu); marker_stop->SetY(sigma); 
            ++point;
        }
        for(int ipu=point; ipu<GetPuType().size(); ++ipu){
            errors[i_seed]->RemovePoint(point);
        }
    }
    TMultiGraph* multi=new TMultiGraph();
    for(auto graph: errors      ){ multi->Add(graph,"lp"); }
    for(auto graph: errors_allPu){ multi->Add(graph,"lp"); }
    multi->Draw("a");
    multi->GetXaxis()->SetTitle(("Mean "+fXTitle).c_str());
    multi->GetYaxis()->SetTitle(("#sigma "+fXTitle).c_str());
    multi->GetYaxis()->SetRangeUser(0,60);
    multi->GetXaxis()->SetRangeUser(50,130);
    leg->Draw();
    seed_leg->Draw();
    gPad->SetGridx();
    DrawCmsStampTurnon(0);
    std::string outName = Form("%s/fits_%s.%s",this->GetOutDir().c_str(),this->GetOutName().c_str(),this->GetOutExtension().c_str());
    can->SaveAs(outName.c_str());
}

void TL1Turnon::DrawCmsStampTurnon(const double & max)
{
    TLatex * latex(new TLatex());
    latex->SetNDC();
    latex->SetTextFont(42);
    if( this->GetSampleName() == "Data" )
        latex->DrawLatex(0.15,0.92,Form("#bf{CMS} #it{Preliminary} %s",this->GetSampleTitle().c_str()));
    else
        latex->DrawLatex(0.15,0.92,Form("#bf{CMS} #it{Simulation Preliminary} %s",this->GetSampleTitle().c_str()));
    latex->SetTextAlign(31);
    latex->DrawLatex(0.92,0.92,Form("%s (13 TeV)",this->GetRun().c_str()));
    //latex->SetTextAlign(32);
    //latex->DrawLatex(0.82,0.25,this->GetAddMark().c_str());

    if(max!=0){
        double min = fXBins.front();
        //double max = fXBins.back();
        TLine * line(new TLine(min,1.,max,1.));
        line->SetLineStyle(7);
        line->DrawClone();
    }
}

TF1* TL1Turnon::fit(TGraphAsymmErrors * eff, double p50){
    std::vector<std::string> fit_functions;
    fit_functions.emplace_back("0.5*(1+TMath::Erf((x-[0])*[1]))");

    // Fit with an exponentially modified Gaussian (EMG)
    // From Wikipedia ( https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution ) the CDF of an EMG is:
    // CDF = \Phi (u,0,v)-e^{-u+v^{2}/2+\log(\Phi (u,v^{2},v))}}
    // \Phi (x,\mu ,\sigma ) is the CDF of a Gaussian distribution,
    // u=\lambda (x-\mu ) 
    // v=\lambda \sigma 
    // \lambda (>0) := the exponential decay parameter
    // \mu := the mean of the Gaussian component
    // \sigma^2 (>0):= the variance of the Gaussian component
    // Which simplifies to:
    // std::string func = "(1+TMath::Erf( (x-[0])*[2]/([1]*[1]))) - exp(-(x - [0] - 0.5/[2]*[1]*[1])/[2])*(1 + TMath::Erf( (x-[0])*[2]/([1]*[1])-1 ))";
    // [0] = \mu, [1] = \sigma, [2] = 1 / \lambda

    // [0] = \mu,  [1] =  1/( \lambda*\sigma^2 ),  [2] = \lambda
    std::string scaled_x  = "(x - [0])*[1]";
    std::string term_1    = Form("0.5 * (1 + TMath::Erf( %s ) )"   , scaled_x.c_str());
    std::string exp_modif = Form("exp(- [2]/[1]*( %s -0.5) )"      , scaled_x.c_str());
    std::string term_2    = Form("0.5 * (1 + TMath::Erf( %s -1) )" , scaled_x.c_str());
    std::string func      = Form("%s - %s*%s",term_1.c_str(),exp_modif.c_str(),term_2.c_str());
    fit_functions.emplace_back(func);


    //std::string func = "[0]*0.5*exp([0]*0.5*(2.0*[1]+[0]*[2]*[2]-2.0*x))*(1-TMath::Erf(([1]+[0]*[2]*[2]-x)/(sqrt(2.0)*[2]))";
    //
    //std::string func = "[0]*0.5*exp([0]*0.5*(2.0*[1]+[0]*[2]*[2]-2.0*x))*(1-TMath::Erf([1]/(sqrt(2)*[2])+[0]*[2]/sqrt(2)-x/(sqrt(2)*[2])))";
    //
    //std::string func = "0.5*(1+[2]*TMath::Gaus(x,[1],1/[0]))*(1+TMath::Erf((x-[1])*[0]))";
    //std::string func = "0.5*(1+TMath::Erf((x-[1])*[0]))";
    //std::string func = "0.5*(1+TMath::Erf((x-[0])*[1])) + [2]*TMath::Gaus( (x-[0])*[1])";
    //std::string func = "0.5*((x<[1])?TMath::Erfc(-(x-[1])*[0]):1+TMath::Erf((x-[1])*[0]))";
    //std::string func = "0.5*((x<[1])?TMath::Erfc(-(x-[1]+log([2]*x))*[0]):1+TMath::Erf((x-[1]+log([2]*x))*[0]))";
    //TF1 fitFcn(Form("fit_%s",eff->GetName()),func.c_str(),fXBins.front(),fXBins.back());
    int count=0;
    std::vector<TF1*> functions;
    for(auto func:fit_functions){
        TF1* fitFcn=new TF1(Form("fit_%s_%d",eff->GetName(),count),func.c_str(),fXBins.front(),fXBins.back());
        functions.emplace_back(fitFcn);

        if(count==0){
            const double mu = p50;
            const double sigma = 10;
            fitFcn->SetParameters(mu,1/sigma);
        }else if(count==1){
            const double mu = functions[0]->GetParameter(0);
            const double sigma = 1/functions[0]->GetParameter(1);
            const double lambda = 0.05; // should be within 0.04 and 0.06 it seems

            const double p0 = mu; 
            const double p1 = 1/(sigma);
            //const double p1 = 1/(sigma * sigma *lambda);
            const double p2 = lambda;
            fitFcn->SetParameters(p0,p1,p2,0);
        }

        //fitFcn.SetParLimits( 2, 0.02,5); // might be an issue since lambda is mixed into p1 and p2

        //fitFcn.SetParameters( 1/150.0,p50 ,0.02);
        //fitFcn.SetParLimits( 2, -1,1);
        //fitFcn.FixParameter(1,p50);
        //fitFcn.SetParameters( 1.000,150.0,(double)p50 );

        TFitResultPtr success=eff->Fit(fitFcn->GetName(),"ESMQ ROB EX0+"); 
        //if((int)success !=0){
        //    cout<<"Fit failed: "<<fitFcn.GetParameter( 0)<<" "<<fitFcn.GetParameter(1)<<endl;
        //    fitFcn.SetParameters( 0,0 );
        //}

        //for(int i=0; i<10; ++i)
        //    //eff->Fit(fitFcn.GetName(),"E0M");
        //    eff->Fit(fitFcn.GetName(),"ELMN"); 

        fitFcn->SetLineColor(eff->GetLineColor());
        fitFcn->SetLineWidth(2);
        fitFcn->SetLineStyle(count);
        TF1* graph_line=dynamic_cast<TF1*>(eff->GetListOfFunctions()->Last());
        if(graph_line){
            graph_line->SetLineColor(eff->GetLineColor());
            graph_line->SetLineWidth(2);
            graph_line->SetLineStyle(2-count);
        }

        ++count;
    }
    eff->GetListOfFunctions()->Print();

    return functions.back();
}

TGraphAsymmErrors GetEfficiency(TH1F * total, TH1F * pass)
{
    TEfficiency * eff = new TEfficiency(*pass, *total);
    std::vector<double> x, y, exl, exh, eyl, eyh;
    double binWidth(0.0);
    for(int bin=1; bin<=total->GetNbinsX(); ++bin)
    {
        binWidth = 0.5*total->GetBinWidth(bin);
        x.push_back(total->GetBinCenter(bin));
        y.push_back(eff->GetEfficiency(bin));
        exl.push_back(binWidth);
        exh.push_back(binWidth);
        eyl.push_back(eff->GetEfficiencyErrorLow(bin));
        eyh.push_back(eff->GetEfficiencyErrorUp(bin));
    }
    x.push_back(total->GetBinCenter(total->GetNbinsX())+2*binWidth);
    y.push_back(eff->GetEfficiency(total->GetNbinsX()+1));
    exl.push_back(binWidth);
    exh.push_back(binWidth);
    eyl.push_back(eff->GetEfficiencyErrorLow(total->GetNbinsX()+1));
    eyh.push_back(eff->GetEfficiencyErrorUp(total->GetNbinsX()+1));

    TGraphAsymmErrors efficiency(x.size(),&(x[0]),&(y[0]),&(exl[0]),&(exh[0]),&(eyl[0]),&(eyh[0]));
    efficiency.SetName(Form("%s_DIV_%s",pass->GetName(),total->GetName()));
    return efficiency;
}

#endif
