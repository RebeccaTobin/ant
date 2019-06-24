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
    const BinSettings mass_bins(2000, 0, 1100);

    //h_TaggerTime = HistFac.makeTH1D("Tagger Time",    // title
    //                                "t [ns]","#",     // xlabel, ylabel
    //                                time_bins, // our binnings
    //                                "h_TaggerTime"    // ROOT object name, auto-generated if omitted
    //                                );
    //h_TaggerCBSubtaction = HistFac.makeTH1D("Tagger CB Subtaction",
    //                                "t [ns]","#",
    //                                time_bins,
    //                                "h_TaggerCBSubtaction"
    //                                );
    h_PromptRandomWithTriggerSimulation = HistFac.makeTH1D("PromptRandom with TriggerSimulation",
                                    "t [ns]","#",
                                    time_bins,
                                    "h_PromptRandomWithTriggerSimulation"
                                    );
    h_MissingMass = HistFac.makeTH1D("Missing Mass",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass"
                                     );
    // Prompt and random windows
    promptrandom.AddPromptRange({ 9, 19 }); // in nanoseconds
    promptrandom.AddRandomRange({ -200, 7 });
    promptrandom.AddRandomRange({ 21, 200 });
}

double Compton::GetMissingMass(const double& incoming_ph_energy,
                            const double& scattered_ph_energy,
                            const double& theta) {
    // Calculates missing mass given the incoming photon energy,
    // scattered photon energy, and scattered photon angle
    return (incoming_ph_energy * scattered_ph_energy)*(1.0 - cos(theta))
            /(incoming_ph_energy - scattered_ph_energy);
}

void Compton::ProcessEvent(const TEvent& event, manager_t&)
{
    // Runs ProcessEvent function in TriggerSimulation file which
    // does the calculations
    triggersimu.ProcessEvent(event);

    const auto target_vec = LorentzVec({0,0,0},ParticleTypeDatabase::Proton.Mass());

    //const auto scattered_ph_vec = LorentzVec({});

    cout << target_vec << endl;

    for (const auto& taggerhit : event.Reconstructed().TaggerHits) {

        // Lorentz momentum vector for incoming photon.
        // Momentum only in the z-direction
        const auto incoming_ph_vec = LorentzVec({0,0,taggerhit.PhotonEnergy},
                                                taggerhit.PhotonEnergy);

        //h_TaggerTime->Fill(taggerhit.Time);

        // Plot tagger time with weighted CBtime subtraction
        // Use this plot to set prompt and random range
        const auto& CorrectedTaggerTime = triggersimu.GetCorrectedTaggerTime(taggerhit);
        //h_TaggerCBSubtaction->Fill(CorrectedTaggerTime);

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
        h_PromptRandomWithTriggerSimulation->Fill(taggerhit.Time, weight);

        // MissingMass plot
        for (const auto& candidatehit : event.Reconstructed().Candidates) {
            double MissingMass = GetMissingMass(taggerhit.PhotonEnergy,
                                         candidatehit.CaloEnergy, candidatehit.Theta);
            h_MissingMass->Fill(MissingMass);


        }
    }
}


void Compton::ShowResult()
{
    ant::canvas(GetName()+": Tagger Time Stuff")
            //<< h_TaggerCBSubtaction
            //<< h_TaggerTime
            << h_PromptRandomWithTriggerSimulation
            << h_MissingMass
            << endc; // actually draws the canvas
}

// A macro that registers the Compton class with Ant so that
// you can call this class using "Ant -p Compton"
AUTO_REGISTER_PHYSICS(Compton)
