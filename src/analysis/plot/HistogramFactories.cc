#include "plot/HistogramFactories.h"

#include "base/std_ext/string.h"

#include "TMath.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace ant;
using namespace ant::analysis;
using namespace ant::std_ext;
using namespace std;


void HistogramFactory::goto_dir() const
{
    if(my_directory)
        my_directory->cd();
}

string HistogramFactory::MakeTitle(const string& title) const
{
    if(title_prefix.empty())
        return title;
    return std_ext::formatter() << title_prefix << ": " << title;
}

TDirectory *HistogramFactory::mkDirNumbered(const string &name, TDirectory *rootdir)
{
    TDirectory* dir = nullptr;

    unsigned n=0;
    do {
        const string dn = (n!=0) ? name+to_string(n) : name;
        ++n;

        rootdir->GetObject(dn.c_str(), dir);

        if(!dir) {

            dir = rootdir->mkdir(dn.c_str());

            if(!dir)
                throw("Can't create output directory \"" + dn +"\"");
        } else {
            dir = nullptr;
        }

    } while(dir==nullptr);

    return dir;
}

string HistogramFactory::GetNextHistName(const string &name) const
{
    if(name.empty()) {
        return formatter() << "hist" << setfill('0') << setw(3) << n_unnamed++;
    } else
        return name;
}

HistogramFactory::HistogramFactory(const string &directory_name, TDirectory* root, const string& title_prefix_):
    title_prefix(title_prefix_)
{

    if(!root)
        root=gDirectory;

    my_directory = mkDirNumbered(directory_name, root);

}


HistogramFactory::HistogramFactory(const string& directory_name, const HistogramFactory& parent, const string& title_prefix_)
  : my_directory(),
    title_prefix
    (
        parent.title_prefix.empty() ?
            title_prefix_ :
            (title_prefix_.empty() ? parent.title_prefix :
                                     std_ext::formatter() << parent.title_prefix << ": " << title_prefix_)
    )
{
    my_directory = parent.my_directory->mkdir(directory_name.c_str());
    if(!my_directory)
        my_directory=gDirectory;
}

void HistogramFactory::SetTitlePrefix(const string& title_prefix_)
{
    title_prefix = title_prefix_;
}

string HistogramFactory::GetTitlePrefix() const
{
    return title_prefix;
}

void HistogramFactory::SetDirDescription(const string &desc)
{
    my_directory->SetTitle(desc.c_str());
}

TH1D *HistogramFactory::makeTH1D(const string &title, const string &xlabel, const string &ylabel, const BinSettings &bins, const string &name) const
{
    auto r = make<TH1D>(GetNextHistName(name).c_str(), MakeTitle(title).c_str(),
                        bins.Bins(), bins.Start(), bins.Stop());
    r->SetXTitle(xlabel.c_str());
    r->SetYTitle(ylabel.c_str());
    return r;
}




TH2D *HistogramFactory::makeTH2D(const string &title,
                                 const string &xlabel,
                                 const string &ylabel,
                                 const BinSettings &xbins,
                                 const BinSettings &ybins,
                                 const string &name) const
{
    auto h = make<TH2D>(GetNextHistName(name).c_str(), MakeTitle(title).c_str(),
                         xbins.Bins(), xbins.Start(), xbins.Stop(),
                         ybins.Bins(), ybins.Start(), ybins.Stop());
    h->SetXTitle(xlabel.c_str());
    h->SetYTitle(ylabel.c_str());
    return h;
}

TH3D *HistogramFactory::makeTH3D(const string &title,
                                 const string &xlabel,
                                 const string &ylabel,
                                 const string &zlabel,
                                 const BinSettings &xbins,
                                 const BinSettings &ybins,
                                 const BinSettings &zbins,
                                 const string &name) const
{
    auto h = make<TH3D>(GetNextHistName(name).c_str(), MakeTitle(title).c_str(),
                       xbins.Bins(), xbins.Start(), xbins.Stop(),
                       ybins.Bins(), ybins.Start(), ybins.Stop(),
                       zbins.Bins(), zbins.Start(), zbins.Stop());
    h->SetXTitle(xlabel.c_str());
    h->SetYTitle(ylabel.c_str());
    h->SetZTitle(zlabel.c_str());
    return h;
}

TTree* HistogramFactory::makeTTree(const string& name) const
{
    return make<TTree>(name.c_str(), MakeTitle(name.c_str()).c_str());
}

HistogramFactory::DirStackPush::DirStackPush(const HistogramFactory& hf): dir(gDirectory)
{
    hf.goto_dir();
}

HistogramFactory::DirStackPush::~DirStackPush()
{
    dir->cd();
}

BinSettings BinSettings::RoundToBinSize(const BinSettings& bins, const double binSize) {

    unsigned nOptBins = unsigned(bins.Length() / binSize);

    const auto factor = unsigned(std::round(double(nOptBins) / double(bins.Bins())));

    if(nOptBins >= binSize) {
         nOptBins /= factor;
    } else {
        nOptBins *= factor;
    }

    const auto length = binSize * nOptBins * factor;

    return BinSettings(nOptBins, interval<double>::CenterWidth(bins.Center(), length));
}

BinSettings BinSettings::Make(const vector<double>& x_values)
{
    if(x_values.size()<2)
        throw runtime_error("Given x_values should contain at least 2 elements");
    vector<double> diffs(x_values.size());
    std::adjacent_difference(x_values.begin(), x_values.end(), diffs.begin());
    if(std::find_if(diffs.begin(), diffs.end(), [] (double v) { return v<=0; }) != diffs.end())
        throw runtime_error("Given x_values not strictly monotonic");
    double mean_distance = std::accumulate(std::next(diffs.begin()), diffs.end(), 0.0) / diffs.size();
    return { static_cast<unsigned>(x_values.size()),
                x_values.front() - mean_distance/2.0,
                x_values.back()  + mean_distance/2.0 };
}
