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
    const BinSettings mass_bins(1000, 0, 1100);

    h_PromptRandomWithTriggerSimulation = HistFac.makeTH1D("PromptRandom with TriggerSimulation",
                                    "t [ns]","#",
                                    time_bins,
                                    "h_PromptRandomWithTriggerSimulation"
                                    );
    h_MissingMass = HistFac.makeTH1D("No Filter",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass"
                                     );
    h_MissingMass1 = HistFac.makeTH1D("With weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass1"
                                     );
    h_MissingMass01 = HistFac.makeTH1D("Veto, no weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass01"
                                     );
    h_MissingMass11 = HistFac.makeTH1D("Veto with weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass11"
                                     );
    h_MissingMass002 = HistFac.makeTH1D("2 particles, no veto or weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass002"
                                     );
    h_MissingMass102 = HistFac.makeTH1D("2 particles, no veto, with weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass102"
                                     );
    h_MissingMass012 = HistFac.makeTH1D("2 particles with veto, no weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass012"
                                     );
    h_MissingMass112 = HistFac.makeTH1D("2 particles with veto and weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass112"
                                     );
    h_MissingMass001 = HistFac.makeTH1D("1 particle, no veto or weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass001"
                                     );
    h_MissingMass101 = HistFac.makeTH1D("1 particle, no veto, with weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass101"
                                     );
    h_MissingMass011 = HistFac.makeTH1D("1 particle with veto, no weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass011"
                                     );
    h_MissingMass111 = HistFac.makeTH1D("1 particle with veto and weights",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass111"
                                     );


    // Prompt and random windows
    promptrandom.AddPromptRange({ 9, 19 }); // in nanoseconds
    promptrandom.AddRandomRange({ -200, 7 });
    promptrandom.AddRandomRange({ 21, 200 });
}


void Compton::ProcessEvent(const TEvent& event, manager_t&)
{

    // Runs ProcessEvent function in TriggerSimulation file which
    // does the calculations
    triggersimu.ProcessEvent(event);

    // Momentum 4 vector for target (i.e. stationary proton)
    const auto target_vec = LorentzVec({0,0,0}, ParticleTypeDatabase::Proton.Mass());

    for (const auto& taggerhit : event.Reconstructed().TaggerHits) {

        // Momentum 4 vector for incoming photon.
        // Note: momentum only in the z-direction
        LorentzVec incoming_ph_vec = LorentzVec({0,0,taggerhit.PhotonEnergy},
                                                taggerhit.PhotonEnergy);
        // Decleration of momentum 4 vector for scattered photon
        // and recoil proton
        LorentzVec scattered_ph_vec;
        LorentzVec recoil_pr_vec;
        double pr_mass;

        // Apply trigger simulation to tagger hits
        // This subtracts a weighted time from the CB (see wiki)
        const auto& CorrectedTaggerTime = triggersimu.GetCorrectedTaggerTime(taggerhit);


        // This assigns weights to the TaggerHits based on which
        // time window they fall into
        promptrandom.SetTaggerTime(CorrectedTaggerTime);

        // When the Tagger hit is in neither the prompt or random
        // window, then skip
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        // Plot taggerhits with weights. Weight of 1.0 in
        // prompt, weight of 0.0 in outside and weight of between
        // 0 and 1 in random
        const double weight = promptrandom.FillWeight();
        h_PromptRandomWithTriggerSimulation->Fill(taggerhit.Time, weight);

        for (const auto& candidate : event.Reconstructed().Candidates) {

            // Momentum 4 vector for scattered photon
            scattered_ph_vec = LorentzVec(vec3(candidate),
                                          candidate.CaloEnergy);
            // Calculating the momentum 4 vector for the possible
            // recoil proton
            recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
            // Calculating the mass of the recoil proton from
            // the momentum vector using .M()
            // Should be 938MeV if there was a Compton
            // event involving these 2 photons
            pr_mass = recoil_pr_vec.M();

            // No filter
            h_MissingMass->Fill(pr_mass);
            // Filter 1: weights
            h_MissingMass1->Fill(pr_mass,weight);

            // Filter 2: check if Candidate is a photon.
            // Arbitrary cut off of 0.2 MeV for how much energy the Veto
            // should get if hit is a photon (theoretically should be 0)
            if (candidate.VetoEnergy < .2)
            {
                h_MissingMass01->Fill(pr_mass);
                h_MissingMass11->Fill(pr_mass,weight);
            }

        }
        // Filter 3: check if there were two particles
        // present in the event
        if (event.Reconstructed().Candidates.size() == 2)
        {
            for (const auto& candidate : event.Reconstructed().Candidates)
            {
                // Doing all this again
                scattered_ph_vec = LorentzVec(vec3(candidate),
                                              candidate.CaloEnergy);
                recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
                pr_mass = recoil_pr_vec.M();

                // 2 particles, no Veto, with and without weights
                h_MissingMass002->Fill(pr_mass);
                h_MissingMass102->Fill(pr_mass,weight);

                if (candidate.VetoEnergy < .2)
                {
                    // 2 particles, with Veto, with and without weights
                    h_MissingMass012->Fill(pr_mass);
                    h_MissingMass112->Fill(pr_mass,weight);
                }
            }
        }

        if (event.Reconstructed().Candidates.size() == 1)
        {
            for (const auto& candidate : event.Reconstructed().Candidates)
            {
                // Doing all this again
                scattered_ph_vec = LorentzVec(vec3(candidate),
                                          candidate.CaloEnergy);
                recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
                pr_mass = recoil_pr_vec.M();

                // 1 particle, no Veto, with and without weights
                h_MissingMass001->Fill(pr_mass);
                h_MissingMass101->Fill(pr_mass,weight);

                if (candidate.VetoEnergy < .2)
                {
                    // 1 particle, with Veto, with and without weights
                    h_MissingMass011->Fill(pr_mass);
                    h_MissingMass111->Fill(pr_mass,weight);
                }
            }
        }
    }
}


void Compton::ShowResult()
{
    ant::canvas(GetName()+": Tagger Time Plots")
            << h_PromptRandomWithTriggerSimulation
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Missing Mass Plots")
            << h_MissingMass
            << h_MissingMass1
            << h_MissingMass01
            << h_MissingMass11
            << h_MissingMass002
            << h_MissingMass102
            << h_MissingMass012
            << h_MissingMass112
            << h_MissingMass001
            << h_MissingMass101
            << h_MissingMass011
            << h_MissingMass111
            << endc;
}

// A macro that registers the Compton class with Ant so that
// you can call this class using "Ant -p Compton"
AUTO_REGISTER_PHYSICS(Compton)
