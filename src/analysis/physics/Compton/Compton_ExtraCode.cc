/*
 * Code to go in final Compton class (but I don't want to run every time)
 * /////////////////////////////////////////////////////////////////////
*/

// In consructor, creating histogram
h_TaggerCBSubtaction = HistFac.makeTH1D("Tagger CB Subtaction",
                                "t [ns]","#",
                                time_bins,
                                "h_TaggerCBSubtaction"
                                );

// In ProcessEvent, plot you would use to get Prompt and
// random time windows (to figure out later)
h_TaggerCBSubtaction->Fill(CorrectedTaggerTime);

/*
 * Code I might use or is helpful for other reasons
 * ////////////////////////////////////////////////
*/



// In constructor, creating a histogram
h_TaggerTime = HistFac.makeTH1D("Tagger Time",    // title
                                "t [ns]","#",     // xlabel, ylabel
                                time_bins, // our binnings
                                "h_TaggerTime"    // ROOT object name, auto-generated if omitted
                                );

// Trying to get missing mass using the formula
// Didn't work, but maybe useful to revisit
double Compton::GetMissingMass(const double& incoming_ph_energy,
                            const double& scattered_ph_energy,
                            const double& theta) {
    // Calculates missing mass given the incoming photon energy,
    // scattered photon energy, and scattered photon angle
    return (incoming_ph_energy * scattered_ph_energy)*(1.0 - cos(theta))
            /(incoming_ph_energy - scattered_ph_energy);
}

// In Candidate for loop. Calculating missing mass
double MissingMass = GetMissingMass(taggerhit.PhotonEnergy,
                               candidate.CaloEnergy, candidate.Theta);
h_MissingMass->Fill(MissingMass);
