// If this file is #include more than once, this command tells
// the compiler to read this file only once
#pragma once

#include "physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class Compton : public Physics {
public:
    // Constructor
    Compton(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

    virtual void ShowResult() override;


private:
    TH1D* h_nClusters;
};

}}}
