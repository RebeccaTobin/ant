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
}

// ?? Don't know why there's no variable after the manager_t& ??
void Compton::ProcessEvent(const TEvent& event, manager_t&)
{
    // Want to inspect what event.Reconstructed() looks like
    for (const auto& recontructed : event.Reconstructed()) {

        cout << recontructed << endl;

}

void Compton::ShowResult()
{

}

// A macro that registers the Compton class with Ant so that
// you can call your class using "Ant -p Compton"
AUTO_REGISTER_PHYSICS(Compton)
