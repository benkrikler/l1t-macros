#include <string>
#include <vector>

#include "../Plotting/tdrstyle.C"
#include "../Event/TL1EventClass.h"
#include "../Utilities/TL1Progress.C"
#include "../Utilities/TL1DateTime.C"
#include "../Plotting/TL1Rates.h"

#include "../Config/ntuple_cfg.h"
#include "../Config/jetRates_cfg.h"

#include "../Debug/DebugHandler.h"

// CHUNK = which chunk of root-files to run over
// NFILES = numbers of root-files
// NJOBS = number of jobs to submit
void makeJetRates(const int & CHUNK, const int & NFILES, const int & NJOBS, const bool & COMBINE)
{
    // Check CHUNK < NJOBS
    DebugHandler::ErrorCheck(CHUNK >= NJOBS, "The CHUNK number exceeds the number of jobs", __FILE__, __LINE__);

    // Set ROOT style
    TStyle * myStyle(new TStyle(TDRStyle()));
    SetMyStyle(55, 0.08, myStyle);

    // Get config objects
    ntuple_cfg * dataset = new ntuple_cfg(GetNtuple_cfg());
    dataset->puType = {};
    dataset->puBins = {};
    std::map< std::string, TL1Rates* > rates = jetRates(dataset);

    // Split the files into CHUNKS:
    // 1 to n; n+1 to 2n; 2n+1 to 3n ... NJOBS*n to NFILES
    // n = nFilesPerJob
    // The final job is typically larger than the others (never smaller)
    std::vector<std::string> inDir;
    std::string outDir( dataset->outDir+"_hadd/RatesJet/" );
    if( !COMBINE )
    {
        outDir = dataset->outDir;
        if( NJOBS > 1 ) outDir += Form("_CHUNK%i",CHUNK);
        outDir += "/RatesJet/";

        int nFilesPerJob( NFILES / NJOBS );
        int finalFile( nFilesPerJob*(1+CHUNK) );
        if( CHUNK == NJOBS-1 ) finalFile = NFILES;

        for(int i=1+(CHUNK*nFilesPerJob); i<=finalFile; ++i)
            inDir.push_back(Form(dataset->inFiles.c_str(),i));
    }
    TL1EventClass * event(new TL1EventClass(inDir));

    // Begin
    for(auto it=rates.begin(); it!=rates.end(); ++it)
    {
        it->second->SetSample(dataset->sampleName, dataset->sampleTitle);
        it->second->SetTrigger(dataset->triggerName, dataset->triggerTitle);
        it->second->SetRun(dataset->run);
        it->second->SetOutDir(outDir);
        if( !COMBINE ) it->second->InitPlots();
        else it->second->OverwritePlots();
    }

    // Loop
    unsigned NEntries(0);
    if( !COMBINE ) NEntries = event->GetPEvent()->GetNEntries();
    while( event->Next() && !COMBINE )
    {
        unsigned position = event->GetPEvent()->GetPosition()+1;
        TL1Progress::PrintProgressBar(position, NEntries);

        // Get the relevant event parameters
        //int pu = event->GetPEvent()->fVertex->nVtx;
        auto l1JetEts = event->fL1JetEt;
        double maxL1JetEt(0.0);
        for(auto it=l1JetEts.begin(); it!=l1JetEts.end(); ++it)
            if( *it > maxL1JetEt )
                maxL1JetEt = *it;
        int nL1Jets = l1JetEts.size();

        if( nL1Jets >= 1 )
            if( rates.find("singleJets") != rates.end() )
                rates["singleJets"]->Fill(maxL1JetEt, 0.0);
        if( nL1Jets >= 2 )
            if( rates.find("doubleJets") != rates.end() )
                rates["doubleJets"]->Fill(maxL1JetEt, 0.0);
        if( nL1Jets >= 3 )
            if( rates.find("tripleJets") != rates.end() )
                rates["tripleJets"]->Fill(maxL1JetEt, 0.0);
        if( nL1Jets >= 4 )
            if( rates.find("quadJets") != rates.end() )
                rates["quadJets"]->Fill(maxL1JetEt, 0.0);
    }

    for(auto it=rates.begin(); it!=rates.end(); ++it)
        it->second->DrawPlots();
    std::cout << "Output saved in:\n\t" << outDir << std::endl;
}