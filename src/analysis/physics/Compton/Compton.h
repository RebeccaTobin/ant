// If this file is #include more than once, this command tells
// the compiler to read this file only once
#pragma once

#include "physics/Physics.h"
// To subtact out random tagger hits
#include "plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"

// A heirarchy of namespaces which generally resembles the
// folder heirarchy
namespace ant {
namespace analysis {
namespace physics {

// Creating a new class called Compton that inherits
// the members of the Physics class. The Physics class is
// defined in "physics.h"

class Compton : public Physics {
public:
    // Constructor created but not defined
    Compton(const std::string& name, OptionsPtr opts);

    // ProcessEvent is a funtion that is used by every physics
    // class. override is used so that when the compiler gets
    // to this code, the ProcessEvent defined here becomes the
    // default. TEvent and manager_t are data types that are
    // defined elswhere in ant (like int or char or double).
    // & indicateds that the variables event and manager are
    // being passed as a reference (which means they will be
    // modified by the function.
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

    virtual void ShowResult() override;


private:
    // TH1 is a root command for making histograms. The D at the
    // end stands for double and it indicates what the height of
    // the bins will be. h_nClusters is the name of the histogram.
    // You would put all your histograms here.
    //TH1D* h_TaggerandClusterTime;
    //TH1D* h_TaggerTime;
    //TH1D* h_TaggerCBSubtaction;
    TH1D* h_PromptRandomWithTriggerSimulation;
    TH1D* h_MissingMass;
    TH1D* h_MissingMass1;
    TH1D* h_MissingMass01;
    TH1D* h_MissingMass11;
    TH1D* h_MissingMass002;
    TH1D* h_MissingMass102;
    TH1D* h_MissingMass012;
    TH1D* h_MissingMass112;
    TH1D* h_MissingMass001;
    TH1D* h_MissingMass101;
    TH1D* h_MissingMass011;
    TH1D* h_MissingMass111;
    PromptRandom::Switch promptrandom;
    utils::TriggerSimulation triggersimu;
};

}}}
