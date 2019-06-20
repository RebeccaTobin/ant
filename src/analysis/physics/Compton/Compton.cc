#include "Compton.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

// Defining the contructor for the Compton class
Compton::Compton(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    // Histograms are created here but not filled

    const BinSettings time_bins(2000, -200, 200);

    h_TaggerTime = HistFac.makeTH1D("Tagger Time",    // title
                                    "t [ns]","#",     // xlabel, ylabel
                                    time_bins, // our binnings
                                    "h_TaggerTime"    // ROOT object name, auto-generated if omitted
                                    );
    h_TaggerCBSubtaction = HistFac.makeTH1D("Tagger CB Subtaction",    // title
                                    "t [ns]","#",     // xlabel, ylabel
                                    time_bins, // our binnings
                                    "h_TaggerCBSubtaction"    // ROOT object name, auto-generated if omitted
                                    );
    h_PromptRandomWithTriggerSimulation = HistFac.makeTH1D("PromptRandom with TriggerSimulation",    // title
                                    "t [ns]","#",     // xlabel, ylabel
                                    time_bins, // our binnings
                                    "h_PromptRandomWithTriggerSimulation"    // ROOT object name, auto-generated if omitted
                                    );
    // Prompt and random windows
    promptrandom.AddPromptRange({ 4,   22}); // in nanoseconds
    promptrandom.AddRandomRange({-200, 0});
    promptrandom.AddRandomRange({ 27,  200});
}

void Compton::ProcessEvent(const TEvent& event, manager_t&)
{
    for (const auto& taggerhit : event.Reconstructed().TaggerHits) {
        h_TaggerTime->Fill(taggerhit.Time);
    }

    // Runs ProcessEvent function in TriggerSimulation file which
    // does the calculations
    triggersimu.ProcessEvent(event);

    for (const auto& clusterhit : event.Reconstructed().Clusters) {
        for (const auto& taggerhit : event.Reconstructed().TaggerHits) {
            h_TaggerCBSubtaction->Fill(taggerhit.Time - clusterhit.Time);
        }
    }

    for (const auto& taggerhit : event.Reconstructed().TaggerHits) {
        // Pass through corrected tagger time into promptrandom

        // Plot tagger time with weighted CBtime subtraction
        // Use this plot to set prompt and random range
        const auto& CorrectedTaggerTime = triggersimu.GetCorrectedTaggerTime(taggerhit);
        h_TaggerCBSubtaction->Fill(CorrectedTaggerTime);

        // This assigns weights to the TaggerHits based on which
        // time window they fall into
        promptrandom.SetTaggerTime(CorrectedTaggerTime);

        // When the Tagger hit is in neither the prompt or random
        // window, then skip
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        // Plot taggerhits with weights. Don't know why taggerhit.Time
        // is used rather than CorrectedTaggerTime
        const double weight = promptrandom.FillWeight();
        h_PromptRandomWithTriggerSimulation->Fill(CorrectedTaggerTime, weight);
    }
}

void Compton::ShowResult()
{
    ant::canvas(GetName()+": Tagger Time Stuff")
            << h_TaggerCBSubtaction
            << h_TaggerTime
            << h_PromptRandomWithTriggerSimulation
            << endc; // actually draws the canvas
}

// A macro that registers the Compton class with Ant so that
// you can call this class using "Ant -p Compton"
AUTO_REGISTER_PHYSICS(Compton)
