#pragma once

#include "analysis/physics/Physics.h"

#include <string>
#include <memory>
#include <vector>


namespace ant {
namespace analysis {
namespace physics {

class EventDisplayTree : public Physics {
protected:
    unsigned n = 0;

    double tapsZ = {};

    static constexpr int clusterMarker = 25;  // Empty square
    static constexpr int trueMarker    = 24;  // Empty circle

    std::vector<data::CandidatePtr> taps_cands;

public:
    EventDisplayTree(const std::string& name, PhysOptPtr opts);
    virtual ~EventDisplayTree();

    void ProcessEvent(const data::Event &event) override;
};

}
}
}
