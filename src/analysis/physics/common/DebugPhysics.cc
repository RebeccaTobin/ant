#include "DebugPhysics.h"

#include "base/Logger.h"
#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

DebugPhysics::DebugPhysics(PhysOptPtr opts): Physics("DebugPhysics", opts) {}

DebugPhysics::~DebugPhysics() {}

void DebugPhysics::ProcessEvent(const data::Event& event)
{
    LOG(INFO) << event;
}

void DebugPhysics::Finish()
{
    LOG(INFO) << "Nop";
}

void DebugPhysics::ShowResult()
{
    LOG(INFO) << "Nop";
}

void DebugPhysics::Initialize(data::Slowcontrol& slowcontrol)
{
    slowcontrol.FaradayCup.Request();
}



DebugPIDAlignment::DebugPIDAlignment(PhysOptPtr opts):
    Physics("DebugPIDAlignment", opts)
{
    const BinSettings bins(360,-180,180);
    angles = HistFac.makeTH2D("PID Phi angles","MCTrue #phi","Rec #phi",bins,bins);
}

DebugPIDAlignment::~DebugPIDAlignment()
{

}

void DebugPIDAlignment::ProcessEvent(const data::Event& event)
{
    if(event.MCTrue().Particles().GetAll().size() == 1) {
        const auto mctrue_phi = event.MCTrue().Particles().GetAll().front()->Phi() * TMath::RadToDeg();

        for(const data::CandidatePtr& c : event.Reconstructed().Candidates()) {
            angles->Fill(mctrue_phi, c->Phi()* TMath::RadToDeg());
        }

        for(const data::Cluster& c : event.Reconstructed().InsaneClusters()) {
            if(c.Detector == Detector_t::Type_t::PID) {
                angles->Fill(mctrue_phi, c.pos.Phi()* TMath::RadToDeg());
            }
        }
    }
}

void DebugPIDAlignment::ShowResult()
{
    canvas("PID angles") << drawoption("colz") << angles << endc;

}

AUTO_REGISTER_PHYSICS(DebugPhysics, "DebugPhysics")
AUTO_REGISTER_PHYSICS(DebugPIDAlignment, "DebugPIDAlignment")
