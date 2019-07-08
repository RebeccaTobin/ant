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
    BinSettings mass_bins(1000, 700, 1100);
    const BinSettings angle_bins1(60 , 120 , 240);
    const BinSettings angle_bins2(90 , 0 , 180);

    h_PromptRandomWithTriggerSimulation = HistFac.makeTH1D("PromptRandom with TriggerSimulation",
                                    "t [ns]","#",
                                    time_bins,
                                    "h_PromptRandomWithTriggerSimulation"
                                    );
    h_MissingMass = HistFac.makeTH1D("All Taggerhits, all Candidates",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass"
                                     );
    h_MissingMass1 = HistFac.makeTH1D("Weighted Taggerhits, all Candidates",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass1"
                                     );
    h_MissingMass01 = HistFac.makeTH1D("All Taggerhits, no Veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass01"
                                     );
    h_MissingMass11 = HistFac.makeTH1D("Weighted Taggerhits, no Veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass11"
                                     );
    h_MissingMass001 = HistFac.makeTH1D("All Taggerhits, 1 particle event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass001"
                                     );
    h_MissingMass101 = HistFac.makeTH1D("Weighted Taggerhits, 1 particle event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass101"
                                     );
    h_MissingMass011 = HistFac.makeTH1D("All Taggerhits, 1 particle event with no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass011"
                                     );
    h_MissingMass111 = HistFac.makeTH1D("Weighted Taggerhits, 1 particle event with no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass111"
                                     );
    h_MissingMass002 = HistFac.makeTH1D("All Taggerhis, 2 particle event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass002"
                                     );
    h_MissingMass102 = HistFac.makeTH1D("Weighted Taggerhis, 2 particle event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass102"
                                     );
    h_MissingMass012 = HistFac.makeTH1D("All Taggerhits, 2 particle event, only one has no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass012"
                                     );
    h_MissingMass112 = HistFac.makeTH1D("Weighted Taggerhits, 2 particle event, only one has no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass112"
                                     );
    h_MissingMass0021 = HistFac.makeTH1D("All Taggerhits, 2 particle event, only plot particle with "
                                         "closer missing mass",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass0021"
                                     );
    h_MissingMass1021 = HistFac.makeTH1D("Weighted Taggerhits, 2 particle event, only plot particle with "
                                         "closer missing mass",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass1021"
                                     );
    h_MissingMass00201 = HistFac.makeTH1D("All Taggerhits, 2 particle coplanar event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass00201"
                                     );
    h_MissingMass10201 = HistFac.makeTH1D("Weighted Taggerhits, 2 particle coplanar event",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass10201"
                                     );
    h_MissingMass01201 = HistFac.makeTH1D("All Taggerhits, 2 particle coplanar event, no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass01201"
                                     );
    h_MissingMass11201 = HistFac.makeTH1D("Weighted Taggerhits, 2 particle coplanar event, no veto",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass11201"
                                     );
    h_MissingMass00211 = HistFac.makeTH1D("All Taggerhits, 2 particle coplanar event,"
                                          " plot only particle with closer missing mass",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass00211"
                                     );
    h_MissingMass10211 = HistFac.makeTH1D("Weighted Taggerhits, 2 particle coplanar event,"
                                          " plot only particle with closer missing mass",
                                     "mass [MeV/c^2]","#",
                                     mass_bins,
                                     "h_MissingMass10211"
                                     );

    h_CoplanarAngle = HistFac.makeTH1D("Coplanar Angle",
                                     "anlge [degrees]","#",
                                     angle_bins1,
                                     "h_CoplanarAngle"
                                     );
    h_OpeningAngle = HistFac.makeTH1D("Opening Angle",
                                     "anlge [degrees]","#",
                                     angle_bins2,
                                     "h_OpeningAngle"
                                     );


    // Prompt and random windows
    promptrandom.AddPromptRange({ 9, 19 }); // in nanoseconds
    promptrandom.AddRandomRange({ -200, 7 });
    promptrandom.AddRandomRange({ 21, 200 });

    tagger_energy_low = 0;
    tagger_energy_high = 1000;
    if (opts->HasOption("low"))
        tagger_energy_low = opts->Get<double>("low", 0);
    if (opts->HasOption("high"))
        tagger_energy_high = opts->Get<double>("high", 1000);

}
// Checks if veto_energy meets threshold for particle
// being considered charged. Returns true if particle
// is charged and false if it is not
bool Compton::IsParticleCharged(double veto_energy)
{
    if (veto_energy < .2) { return false; }
    else { return true; }
}

// Checks if two particles are a photon and proton
// based on their veto energy. Returns 0 is they are
// not, returns the energy of the photon if they are
int Compton::IsPhotonProton(TCandidateList candidates)
{
    bool is1charged = true;
    bool is2charged = true;

    if (candidates.front().VetoEnergy < .2)
    {
        is1charged = false;
    }

    if (candidates.back().VetoEnergy < .2)
    {
        is2charged = false;
    }

    if ((is1charged == false) & (is2charged == true))
    {
        return candidates.front();
    }
    else if ((is1charged == true) & (is2charged == false))
    {
        return candidates.back();
    }

    else { return 0; }
}

// Input: a candidate and the 4 momentum vectors the the
// incoming photon and proton target. Output: the missing
// mass
double GetMissingMass(const TCandidate candidate,
                 LorentzVec target_vec, LorentzVec incoming_ph_vec)
{
    LorentzVec scattered_ph_vec;
    LorentzVec recoil_pr_vec;

    // Momentum 4 vector for scattered photon
    scattered_ph_vec = LorentzVec(vec3(candidate),
                                  candidate.CaloEnergy);
    // Calculating the momentum 4 vector for the possible
    // recoil proton
    recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
    // Calculating the mass of the recoil proton from
    // the 4 momentum vector using .M()
    // Should be 938MeV if there was a Compton
    // event involving these 2 photons
    return recoil_pr_vec.M();
}

void Compton::ProcessEvent(const TEvent& event, manager_t&)
{

    // Runs ProcessEvent function in TriggerSimulation file which
    // does the calculations
    triggersimu.ProcessEvent(event);

    // Momentum 4 vector for target (i.e. stationary proton)
    const auto target_vec = LorentzVec({0,0,0}, ParticleTypeDatabase::Proton.Mass());

    for (const auto& taggerhit : event.Reconstructed().TaggerHits)
    {
        if ((taggerhit.PhotonEnergy < tagger_energy_low) ||
                (taggerhit.PhotonEnergy > tagger_energy_high))
        {
            continue;
        }


        // Momentum 4 vector for incoming photon.
        // Note: momentum only in the z-direction
        LorentzVec incoming_ph_vec = LorentzVec({0,0,taggerhit.PhotonEnergy},
                                                taggerhit.PhotonEnergy);
        // Decleration of momentum 4 vector for scattered photon
        // and recoil proton
        LorentzVec scattered_ph_vec;
        LorentzVec recoil_pr_vec;
        double missing_mass;

        // Apply trigger simulation to tagger hits
        // This subtracts a weighted time from the CB (see wiki)
        const auto& CorrectedTaggerTime =
                triggersimu.GetCorrectedTaggerTime(taggerhit);


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

            missing_mass = GetMissingMass(candidate, target_vec,incoming_ph_vec);

            h_MissingMass->Fill(missing_mass);
            h_MissingMass1->Fill(missing_mass,weight);

            // Filter 2: check if Candidate is uncharged (could it be a photon?)
            if (Compton::IsParticleCharged(candidate.VetoEnergy) == false)
            {
                h_MissingMass01->Fill(pr_mass);
                h_MissingMass11->Fill(pr_mass,weight);
            }

            // Opening angle
            double open_ang = recoil_pr_vec.Angle(scattered_ph_vec);
            h_OpeningAngle->Fill(180*open_ang/M_PI);


        }
        // Filter 3: check if there were two particles
        // present in the event
        if (event.Reconstructed().Candidates.size() == 2)
        {
            // Using both candidates to calc and plot missing mass
            for (const auto& candidate : event.Reconstructed().Candidates)
            {
                // Doing all this again
                scattered_ph_vec = LorentzVec(vec3(candidate),
                                              candidate.CaloEnergy);
                recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
                pr_mass = recoil_pr_vec.M();

                // 2 particles in event, with and without weights
                h_MissingMass002->Fill(pr_mass);
                h_MissingMass102->Fill(pr_mass,weight);

            }

            photon_energy = IsPhotonProton(event.Reconstructed().Candidates)

            if (photon_energy != 0)
            {
                missing_mass = GetMissingMass(Candidates.front(),
                                              target_vec,incoming_ph_vec);
                h_MissingMass012->Fill(missing_mass);
                h_MissingMass112->Fill(missing_mass,weight);

            }
            if (IsPhotonProton(event.Reconstructed().Candidates.front().VetoEnergy
                               ,event.Reconstructed().Candidates.back().VetoEnergy) == 2)
            {
                // Doing all this again
                scattered_ph_vec = LorentzVec(vec3(event.Reconstructed().
                                              Candidates.back()),
                                              event.Reconstructed().Candidates.
                                              back().CaloEnergy);
                recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
                pr_mass = recoil_pr_vec.M();
                h_MissingMass012->Fill(pr_mass);
                h_MissingMass112->Fill(pr_mass,weight);
            }

            // Closer missing mass
            double front_energy =
                    event.Reconstructed().Candidates.front().CaloEnergy;
            double back_energy =
                    event.Reconstructed().Candidates.back().CaloEnergy;

            if (abs(front_energy - 938) < abs(back_energy - 938))
            {
                // Doing all this again
                scattered_ph_vec = LorentzVec(vec3(event.Reconstructed().
                                                   Candidates.front()),
                                                   front_energy);
                recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
                pr_mass = recoil_pr_vec.M();
                h_MissingMass0021->Fill(pr_mass);
                h_MissingMass1021->Fill(pr_mass,weight);
            }
            else if (abs(front_energy - 938) > abs(back_energy - 938))
            {
                // Doing all this again
                scattered_ph_vec = LorentzVec(vec3(event.Reconstructed().
                                                   Candidates.back()),
                                                   back_energy);
                recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
                pr_mass = recoil_pr_vec.M();
                h_MissingMass0021->Fill(pr_mass);
                h_MissingMass1021->Fill(pr_mass,weight);
            }
            else { continue; }

            // Coplanar filter
            double front_phi = event.Reconstructed().Candidates.front().Phi;
            double back_phi = event.Reconstructed().Candidates.back().Phi;

            double diff;

            if (front_phi > back_phi)
            {
                diff = front_phi - back_phi;
            }
            if (back_phi > front_phi)
            {
                diff = back_phi - front_phi;
            }

            // plus or minus 10 degrees from 180
            if (((180*diff/M_PI) > 170) & ((180*diff/M_PI) < 190))
            {
                scattered_ph_vec = LorentzVec(vec3(event.Reconstructed().
                                                   Candidates.front()),
                                                   front_energy);
                recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
                pr_mass = recoil_pr_vec.M();
                h_MissingMass00201->Fill(pr_mass);
                h_MissingMass10201->Fill(pr_mass,weight);

                // Veto filter
                if (IsPhotonProton(event.Reconstructed().Candidates.front().VetoEnergy
                                   ,event.Reconstructed().Candidates.back().VetoEnergy) == 1)
                {
                    // Doing all this again
                    scattered_ph_vec = LorentzVec(vec3(event.Reconstructed().
                                                  Candidates.front()),
                                                  event.Reconstructed().Candidates.
                                                  front().CaloEnergy);
                    recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
                    pr_mass = recoil_pr_vec.M();
                    h_MissingMass01201->Fill(pr_mass);
                    h_MissingMass11201->Fill(pr_mass,weight);

                }
                if (IsPhotonProton(event.Reconstructed().Candidates.front().VetoEnergy
                                   ,event.Reconstructed().Candidates.back().VetoEnergy) == 2)
                {
                    // Doing all this again
                    scattered_ph_vec = LorentzVec(vec3(event.Reconstructed().
                                                  Candidates.back()),
                                                  event.Reconstructed().Candidates.
                                                  back().CaloEnergy);
                    recoil_pr_vec = incoming_ph_vec + target_vec - scattered_ph_vec;
                    pr_mass = recoil_pr_vec.M();
                    h_MissingMass01201->Fill(pr_mass);
                    h_MissingMass11201->Fill(pr_mass,weight);

                    // Closer missing mass filter

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

                if (Compton::IsParticleCharged(candidate.VetoEnergy) == false)
                {
                    // 1 uncharged particle with and without weights
                    h_MissingMass011->Fill(pr_mass);
                    h_MissingMass111->Fill(pr_mass,weight);
                }
            }
        }
    }
}


void Compton::ShowResult()
{
    ant::canvas(GetName()+": Tagger Time Plots, incoming photon energy range: "
                          )
            << h_PromptRandomWithTriggerSimulation
            << endc; // actually draws the canvas

    ant::canvas(GetName()+": Missing Mass Plots, incoming photon energy range: "
                          )
            << h_MissingMass
            << h_MissingMass1
            << h_MissingMass01
            << h_MissingMass11
            << endc;

    ant::canvas(GetName()+": Missing Mass Plots 2, incoming photon energy range: "
                          )
            << h_MissingMass001
            << h_MissingMass101
            << h_MissingMass011
            << h_MissingMass111
            << h_MissingMass002
            << h_MissingMass102
            << h_MissingMass012
            << h_MissingMass112
            << h_MissingMass0021
            << h_MissingMass1021
            << endc;

    ant::canvas(GetName()+": Coplanar Missing Mass Plots, incoming photon energy range: "
                          )
            << h_MissingMass00201
            << h_MissingMass10201
            << h_MissingMass01201
            << h_MissingMass11201
            << h_MissingMass00211
            << h_MissingMass10211
            << endc;

    ant::canvas(GetName()+": Coplaner and Opening Angle, incoming photon energy range: "
                          )
            << h_CoplanarAngle
            << h_OpeningAngle
            << endc;
}

// A macro that registers the Compton class with Ant so that
// you can call this class using "Ant -p Compton"
AUTO_REGISTER_PHYSICS(Compton)
