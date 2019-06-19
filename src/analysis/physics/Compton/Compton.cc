#include "Compton.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

// Defining the contructor from the Compton class. The base
// class Physics is included with the same parameters
// ?? Not too sure about this formatting ??
Compton::Compton(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    // This stuff will happen automatically when an object
    // in the Compton class is created. "In the constructor
    // we do not do anything so far, we just pass the name
    // and the options to the base class to handle some stuff
    // for us, like setting the name of our class and check
    // for some options. The way the information is passed
    // to the base class with calling the constructor after
    // the colon of our own constructor is called member
    // initializer list. It's just a faster way of initializing
    // the base class compared to something like
    // Tutorial::Turorial(const string& name, OptionsPtr opts) { Physics(name, opts); }
    // and is used nearly everywhere within Ant.""
    // ^don't know what that means

    const BinSettings time_bins(2000, -200, 200);

    //h_TaggerandClusterTime = HistFac.makeTH1D("Subtracted Time",    // title
    //                                "t [ns]","#",     // xlabel, ylabel
    //                                time_bins, // our binnings
    //                                "h_TaggerTime"    // ROOT object name, auto-generated if omitted
    //                                );
    h_TaggerTime = HistFac.makeTH1D("Tagger Time",    // title
                                    "t [ns]","#",     // xlabel, ylabel
                                    time_bins, // our binnings
                                    "h_TaggerTime"    // ROOT object name, auto-generated if omitted
                                    );
    // PromptRandom stuff
    promptrandom.AddPromptRange({ -7,   7}); // in nanoseconds
    promptrandom.AddRandomRange({-50, -10});
    promptrandom.AddRandomRange({ 10,  50});
}

// ?? Don't know why there's no variable after the manager_t& ??
// See wiki for more
void Compton::ProcessEvent(const TEvent& event, manager_t&)
{
    // Using corrected tagger time
    triggersimu.ProcessEvent(event);

    for (const auto& taggerhit : event.Reconstructed().TaggerHits) {
        // Pass through corrected tagger time into promptrandom
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));

        // The Tagger hit is in neither the prompt or random window
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        // A weight it given to the hit depending on whether the hit
        // is in the prompt or random window. This allows the random
        // hits to automatically to subtracted out.
        const double weight = promptrandom.FillWeight();
        h_TaggerTime->Fill(taggerhit.Time, weight);
    }

    // This line loops over all the clusters and taggerhits and
    // subtracts the times at which they occur and plots that
    // subtracted time. The peak region represents combinations
    // of taggerhit times and cluster times in which those events
    // could possibly be correlated.

    //for (const auto& clusterhit : event.Reconstructed().Clusters) {
    //    for (const auto& taggerhit : event.Reconstructed().TaggerHits) {
    //        h_TaggerandClusterTime->Fill(taggerhit.Time - clusterhit.Time);

//        }
//    }
}

void Compton::ShowResult()
{
    ant::canvas(GetName()+": Tagger Time with Prompt Random")
            //<< h_TaggerandClusterTime
            << h_TaggerTime
            << endc; // actually draws the canvas
}

// A macro that registers the Compton class with Ant so that
// you can call your class using "Ant -p Compton"
AUTO_REGISTER_PHYSICS(Compton)
