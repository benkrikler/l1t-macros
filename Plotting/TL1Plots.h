#ifndef TL1PLOTS_H
#define TL1PLOTS_H

#include <string>
#include <stdlib.h>

#include <TH1F.h>
#include <TGraph.h>
#include <TRandom3.h>

#include "../Debug/DebugHandler.h"

class TL1Plots
{
    public:
        TL1Plots();
        ~TL1Plots();

        virtual void InitPlots() = 0;
        virtual void OverwritePlots() = 0;
        virtual void Fill(const double & xVal, const double & yVal, const int & pu) = 0;
        virtual void DrawPlots(const char* name_append=NULL) = 0;
        virtual void NormaliseArea(double norm){ }

        virtual void SetOverwriteNames(const std::string & owRootName, const std::string & owHistName);
        virtual void SetSample(const std::string & sampleName, const std::string & sampleTitle);
        virtual void SetTrigger(const std::string & triggerName, const std::string & triggerTitle);
        virtual void SetRun(const std::string & run);
        virtual void SetOutName(const std::string & outName);
        virtual void SetOutDir(const std::string & outDir);
        virtual void SetOutExtension(const std::string & outExt);
        virtual void SetAddMark(const std::string & addMark);
        virtual void SetPuType(const std::vector<std::string> & puType);
        virtual void SetPuBins(const std::vector<int> & puBins);
        virtual void SetPuFile(const std::string & puFileName);
        void SetColor(TH1 * obj, int pos, int max, bool setFill=false);
        void SetColor(TGraph * obj, int pos, int max, bool setFill=false);
        int CalculateColor(int pos,int max);

        double GetPuWeight(int pu);

        std::string GetDrawOption()const{return fDrawOption;}
        void SetDrawOption(const std::string& opts){ fDrawOption=opts;}

    public:
        std::string GetOverwriteRootFilename() const;
        std::string GetOverwriteHistname() const;
        std::string GetSampleName() const;
        std::string GetTriggerName() const;
        std::string GetRun() const;
        std::string GetSampleTitle() const;
        std::string GetTriggerTitle() const;
        std::string GetOutName() const;
        std::string GetOutDir() const;
        std::string GetOutExtension() const;
        std::string GetAddMark() const;
        std::vector<std::string> GetPuType() const;
        std::vector<int> GetPuBins() const;

        double GetRnd() const;

    private:
        TH1F * fPuWeights;
        TRandom3 * fRnd;

        std::string fOverwriteRootFilename, fOverwriteHistname;
        std::string fSampleName, fTriggerName, fRun;
        std::string fSampleTitle, fTriggerTitle;
        std::string fOutName, fOutDir, fOutExtension;
        std::string fAddMark;
        std::vector<std::string> fPuType;
        std::vector<int> fPuBins;
        std::string fDrawOption;

};

TL1Plots::TL1Plots() :
    fRnd(new TRandom3()),
	fOutExtension("png")
{
}

TL1Plots::~TL1Plots()
{
    delete fPuWeights;
    delete fRnd;
}

void TL1Plots::SetOverwriteNames(const std::string & owRootName, const std::string & owHistName)
{
    fOverwriteRootFilename = owRootName; 
    fOverwriteHistname = owHistName;
}

void TL1Plots::SetSample(const std::string & sampleName, const std::string & sampleTitle)
{
    fSampleName = sampleName;
    fSampleTitle = sampleTitle;
}

void TL1Plots::SetTrigger(const std::string & triggerName, const std::string & triggerTitle)
{
    fTriggerName = triggerName;
    fTriggerTitle = triggerTitle;
}

void TL1Plots::SetRun(const std::string & run)
{
    fRun = run;
}

void TL1Plots::SetOutName(const std::string & outName)
{
    fOutName = outName;
}

void TL1Plots::SetOutDir(const std::string & outDir)
{
    std::system(Form("mkdir -p \"%s\"",outDir.c_str()));
    fOutDir = outDir;
}

void TL1Plots::SetOutExtension(const std::string & outExtension)
{
    fOutExtension = outExtension;
}

void TL1Plots::SetAddMark(const std::string & addMark)
{
    fAddMark = addMark;
}

void TL1Plots::SetPuType(const std::vector<std::string> & puType)
{
    fPuType = puType;
}

void TL1Plots::SetPuBins(const std::vector<int> & puBins)
{
    fPuBins = puBins;
}

void TL1Plots::SetPuFile(const std::string & puFileName)
{
    TFile * fPuFile = TFile::Open(puFileName.c_str(),"READ");
    DebugHandler::CheckTFile(fPuFile, __FILE__, __LINE__);

    fPuWeights = (TH1F*)fPuFile->Get("puRatio");
    fPuWeights->SetDirectory(0);
    delete fPuFile;
}

int TL1Plots::CalculateColor(int pos,int max){
    double modifier(0.15), colorIndex;
    int colour(1);
    double fraction = (double)(pos)/(double)(max-1);

    if( pos > max-1 || pos < 0 || max < 0 ) colour = 1;
    else
    {
        colorIndex = (fraction * (1.0-2.0*modifier) + modifier) * gStyle->GetNumberOfColors();
        colour = gStyle->GetColorPalette(colorIndex);
    }
    return colour;
}

void TL1Plots::SetColor(TH1 * obj, int pos, int max, bool setFill)
{
    int colour=CalculateColor(pos,max);
    obj->SetLineColor(colour);
    obj->SetMarkerColor(colour);
    if(setFill) obj->SetFillColor(colour);
    else        obj->SetFillColor(0);
}

void TL1Plots::SetColor(TGraph * obj, int pos, int max, bool setFill)
{
    int colour=CalculateColor(pos,max);
    obj->SetLineColor(colour);
    obj->SetMarkerColor(colour);
    if(setFill) obj->SetFillColor(colour);
    else        obj->SetFillColor(0);
}

std::string TL1Plots::GetOverwriteRootFilename() const
{
    return fOverwriteRootFilename;
}

std::string TL1Plots::GetOverwriteHistname() const
{
    return fOverwriteHistname;
}

double TL1Plots::GetPuWeight(int pu)
{
    if( this->GetSampleName() == "Data" || pu <= 0 ) return 1.0;
    int bin(fPuWeights->GetXaxis()->FindFixBin(pu));
    return fPuWeights->GetBinContent(bin);
}

std::string TL1Plots::GetSampleName() const
{
    return fSampleName;
}

std::string TL1Plots::GetTriggerName() const
{
    return fTriggerName;
}

std::string TL1Plots::GetRun() const
{
    return fRun;
}

std::string TL1Plots::GetSampleTitle() const
{
    return fSampleTitle;
}

std::string TL1Plots::GetTriggerTitle() const
{
    return fTriggerTitle;
}

std::string TL1Plots::GetOutName() const
{
    return fOutName;
}

std::string TL1Plots::GetOutDir() const
{
    return fOutDir;
}

std::string TL1Plots::GetOutExtension() const
{
    return fOutExtension;
}

std::string TL1Plots::GetAddMark() const
{
    return fAddMark;
}

std::vector<std::string> TL1Plots::GetPuType() const
{
    return fPuType;
}

std::vector<int> TL1Plots::GetPuBins() const
{
    return fPuBins;
}

double TL1Plots::GetRnd() const
{
    return fRnd->Rndm();
}

#endif
