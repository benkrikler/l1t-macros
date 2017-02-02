#include <string>
#include <vector>

#include "../Plotting/tdrstyle.C"
#include "../Event/TL1EventClass.h"
#include "../Utilities/TL1Progress.C"

#include "../Config/ntuple_cfg.h"
#include "../Config/sumRates_cfg.h"

#include "../Debug/DebugHandler.h"

// CHUNK = which chunk of root-files to run over
// NJOBS = number of jobs to submit
void makePlots(const int & CHUNK, const int & NJOBS, const int & NENT, const bool & COMBINE)
{
    // Check CHUNK < NJOBS
    DebugHandler::ErrorCheck(CHUNK >= NJOBS, "The CHUNK number exceeds the number of jobs", __FILE__, __LINE__);

    // Set ROOT style
    TStyle * myStyle(new TStyle(TDRStyle()));
    SetMyStyle(55, 0.08, myStyle);

    // Get config objects
    ntuple_cfg * dataset = new ntuple_cfg(GetNtuple_cfg());
    auto plots = MakePlots(dataset);

    std::vector<std::string> inDir = dataset->inFiles;
    std::string outDir( dataset->outDir+"_hadd/Rates/" );
    if(!COMBINE) outDir = dataset->outDir + Form("_CHUNK%i/Rates/",CHUNK);
    else inDir.clear();
    TL1EventClass * event(new TL1EventClass(inDir));

    // Begin
    for(auto plot: plots)
    {
        plot.second->SetSample(dataset->sampleName, dataset->sampleTitle);
        plot.second->SetTrigger(dataset->triggerName, dataset->triggerTitle);
        plot.second->SetRun(dataset->run);
        plot.second->SetPuType(dataset->puType);
        plot.second->SetPuBins(dataset->puBins);
        plot.second->SetOutDir(outDir);
        if( !COMBINE ) plot.second->InitPlots();
        else plot.second->OverwritePlots();
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
	FillPlots(i,event);
    }

    for(auto plot: plots)
        plot.second->DrawPlots();
    std::cout << "Output saved in:\n\t" << outDir << std::endl;
}

#define MAIN_FUNCTION(macro) \
void macro(const int & CHUNK, const int & NJOBS, const int & NENT, const bool & COMBINE){ \
    makePlots(CHUNK,NJOBS, NENT, COMBINE); \
}
