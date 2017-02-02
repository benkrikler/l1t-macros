#include "../Plotting/TL1RootHists.h"
#include "../Plotting/TL1Rates.h"
#include "../Plotting/TL1Turnon.h"
#include "../Event/TL1EventClass.h"

class ntuple_cfg;
typedef std::map<std::string,TL1Plots*> PlotList;
PlotList MakePlots(ntuple_cfg*);
void FillPlots(int i_entry, const TL1EventClass*);

PlotList MakePlots(ntuple_cfg*){
    PlotList plots;

    // Tower 28 MET
    //plots.

    // MET in all towers

    // MET in all towers except 28

    return plots;
}

void FillPlots(int i_entry, const TL1EventClass*){
}

#include "../MakePlots/makePlots.cxx"
MAIN_FUNCTION(studyTower28MET)
