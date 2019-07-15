#include "GetPromptRandomWindows.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

// Defining the contructor for the Compton class
GetPromptRandomWindows::GetPromptRandomWindows(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    // Histograms are created here but not filled

    const BinSettings time_bins(1000, -100, 100);

    h_TaggerTime = HistFac.makeTH1D("Tagger Time",    // title
                                    "t [ns]","#",     // xlabel, ylabel
                                    time_bins,        // binnings
                                    "h_TaggerTime"    // ROOT object name, auto-generated if omitted
                                    );
}

void GetPromptRandomWindows::ProcessEvent(const TEvent& event, manager_t&)
{
    // Runs ProcessEvent function in TriggerSimulation file which
    // does the calculations
    triggersimu.ProcessEvent(event);

    for (const auto& taggerhit : event.Reconstructed().TaggerHits)
    {

        // Apply trigger simulation to tagger hits
        // This subtracts a weighted time from the CB (see wiki)
        const auto& CorrectedTaggerTime =
                triggersimu.GetCorrectedTaggerTime(taggerhit);

        promptrandom.SetTaggerTime(CorrectedTaggerTime);
        h_TaggerTime->Fill(CorrectedTaggerTime);
    }
}

void GetPromptRandomWindows::ShowResult()
{
    ant::canvas(GetName()+": Tagger Time Plot")
            << h_TaggerTime
            << endc; // actually draws the canvas
}

// A macro that registers the Compton class with Ant so that
// you can call this class using "Ant -p Compton"
AUTO_REGISTER_PHYSICS(GetPromptRandomWindows)
