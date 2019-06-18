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

    const BinSettings tagger_time_bins(2000, -200, 200);

    h_TaggerTime = HistFac.makeTH1D("Tagger Time",    // title
                                    "t [ns]","#",     // xlabel, ylabel
                                    tagger_time_bins, // our binnings
                                    "h_TaggerTime"    // ROOT object name, auto-generated if omitted
                                    );

    h_nClusters = HistFac.makeTH1D("Number of Cluster", "nClusters",
                                   "#", BinSettings(15,0,15), "h_nClusters");
}

// ?? Don't know why there's no variable after the manager_t& ??
// See wiki for more
void Compton::ProcessEvent(const TEvent& event, manager_t&)
{

    // This line loops over the vector TaggerHits (see long form of code
    // below) with parameter taggerhit.
    // for (const auto& taggerhit = event.Reconstructed().TaggerHits.begin();
    // taggerhit != event.Reconstructed().TaggerHits.end(); ++taggerhit)
    for (const auto& taggerhit : event.Reconstructed().TaggerHits) {

        h_TaggerTime->Fill(taggerhit.Time);
    }

    // The number of clusters is independent from the Tagger hits
    // so this line is outside the for loop
    h_nClusters->Fill(event.Reconstructed().Clusters.size());

}

void Compton::ShowResult()
{
    ant::canvas(GetName()+": Basic plots")
            << h_TaggerTime
            << h_nClusters
            << endc; //actually draws the canvas
}

// A macro that registers the Compton class with Ant so that
// you can call your class using "Ant -p Compton"
AUTO_REGISTER_PHYSICS(Compton)
