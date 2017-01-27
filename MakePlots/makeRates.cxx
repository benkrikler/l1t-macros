#include <string>
#include <vector>

#include "../Plotting/tdrstyle.C"
#include "../Event/TL1EventClass.h"
#include "../Utilities/TL1Progress.C"
#include "../Plotting/TL1Rates.h"

#include "../Config/ntuple_cfg.h"
#include "../Config/sumRates_cfg.h"

#include "../Debug/DebugHandler.h"

// CHUNK = which chunk of root-files to run over
// NJOBS = number of jobs to submit
void makeRates(const int & CHUNK, const int & NJOBS, const int & NENT, const bool & COMBINE)
{
    // Check CHUNK < NJOBS
    DebugHandler::ErrorCheck(CHUNK >= NJOBS, "The CHUNK number exceeds the number of jobs", __FILE__, __LINE__);

    // Set ROOT style
    TStyle * myStyle(new TStyle(TDRStyle()));
    SetMyStyle(55, 0.08, myStyle);

    // Get config objects
    ntuple_cfg * dataset = new ntuple_cfg(GetNtuple_cfg());
    std::map< std::string, TL1Rates* > rates = sumRates(dataset);

    std::vector<std::string> inDir = dataset->inFiles;
    std::string outDir( dataset->outDir+"_hadd/Rates/" );
    if(!COMBINE) outDir = dataset->outDir + Form("_CHUNK%i/Rates/",CHUNK);
    else inDir.clear();
    TL1EventClass * event(new TL1EventClass(inDir));

    // Begin
    for(auto it=rates.begin(); it!=rates.end(); ++it)
    {
        it->second->SetSample(dataset->sampleName, dataset->sampleTitle);
        it->second->SetTrigger(dataset->triggerName, dataset->triggerTitle);
        it->second->SetRun(dataset->run);
        it->second->SetPuType(dataset->puType);
        it->second->SetPuBins(dataset->puBins);
        it->second->SetOutDir(outDir);
        if( !COMBINE ) it->second->InitPlots();
        else it->second->OverwritePlots();
    }

    unsigned start(0), end(0);
    if( !COMBINE )
    {
        unsigned NEvents(NENT / NJOBS);
        start = CHUNK * NEvents;
        end   = (CHUNK+1) * NEvents;
        if( CHUNK == NJOBS-1 ) end = NENT;
    }

    // Loop
    for(int i=start; i<end && !COMBINE; ++i)
    {
        event->GetEntry(i);
        TL1Progress::PrintProgressBar(i-start, end-start);

        const int pu = event->GetPEvent()->fVertex->nVtx;
	//if(pu<40) continue;

        // Get the relevant event parameters
        double l1MetBE = event->fL1Met;
	double l1MetBEEmu = event->fL1EmuMet;
	double l1MetBERecalcEmu = event->fRecalcL1EmuMet;
	double l1MetBERecalc = event->fRecalcL1Met;
        double l1MetHF = event->fL1MetHF;

        if( rates.find("l1MetBE") != rates.end() )
            rates["l1MetBE"]->Fill(l1MetBE,1, pu);
	if( rates.find("l1MetBEEmu") != rates.end() )
            rates["l1MetBEEmu"]->Fill(l1MetBEEmu,1, pu);
	if( rates.find("l1MetBERecalc") != rates.end() )
            rates["l1MetBERecalc"]->Fill(l1MetBERecalc,1,pu);
	if( rates.find("l1MetBERecalcEmu") != rates.end() )
	  rates["l1MetBERecalcEmu"]->Fill(l1MetBERecalcEmu,1, pu);
        if( rates.find("l1MetHF") != rates.end() )
            rates["l1MetHF"]->Fill(l1MetHF,1, pu);
    }

    for(auto it=rates.begin(); it!=rates.end(); ++it)
        it->second->DrawPlots();
    std::cout << "Output saved in:\n\t" << outDir << std::endl;
}
