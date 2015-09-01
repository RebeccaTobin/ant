#include "DataBase.h"

#include "base/WrapTFile.h"
#include "base/interval.h"
#include "base/Logger.h"

#include "base/std_ext/misc.h"

#include "TTree.h"

#include <dirent.h>

using namespace std;
using namespace ant;
using namespace ant::calibration;

DataBase::DataBase():
    cm_treename_prefix("calibration-"),
    cm_branchname("cdata")
{

}

bool DataBase::ReadFromFile(const std::string& filename)
{
    try {
        WrapTFileInput dataFile;
        dataFile.OpenFile(filename);

        for( TTree* calibtree: dataFile.GetListOf<TTree>())
        {
            string tname(calibtree->GetName());
            if(!std_ext::begins_with(tname, cm_treename_prefix))
                continue;
            const TCalibrationData* cdata = nullptr;
            calibtree->SetBranchAddress(cm_branchname.c_str(),&cdata);
            for (Long64_t entry = 0; entry < calibtree->GetEntries(); ++entry)
            {
                calibtree->GetEntry(entry);
                DataMap[cdata->CalibrationID].push_back(*cdata); //emplace???
            }
        }
        return true;
    } catch (...) {
        LOG(WARNING) << "Cannot open " << filename;
        return false;
    }

}

bool DataBase::ReadFromFolder(const string& folder)
{
    // open all root files in that folder
    DIR* dir = opendir (folder.c_str());
    if(dir==nullptr)
        return false;
    dirent* ent;
    while((ent = readdir(dir))) {
        string filename = ent->d_name;
        if(!std_ext::string_ends_with(filename,".root"))
            continue;
        ReadFromFile(folder+"/"+filename);
    }
    closedir(dir);
    return true;
}

void DataBase::WriteToTree(WrapTFileOutput& file, const DataMap_t::value_type& calibration) const
{
    string tname = cm_treename_prefix + calibration.first;

    TTree* currentTree = file.CreateInside<TTree>(tname.c_str(),tname.c_str());
    const TCalibrationData* cdataptr = nullptr;
    currentTree->Branch(cm_branchname.c_str(), std::addressof(cdataptr));

    for(const auto& cdata : calibration.second) {
        cdataptr = std::addressof(cdata);
        currentTree->Fill();
    }
}

void DataBase::WriteToFile(const std::string& filename) const
{
    // write everything into one single file
    WrapTFileOutput file(filename,
                         WrapTFileOutput::mode_t::recreate);
    // loop over map and write a new tree for each calibrationID
    for(const auto& calibration: DataMap)
    {
        WriteToTree(file, calibration);
    }
}

void DataBase::WriteToFolder(const string& folder) const
{
    // we cannot use a simple Update mode of the file
    // since we want to recreate the files from the database
    // (items might have been deleted...)

    map<string, shared_ptr<WrapTFileOutput>> files_by_filename;
    map<string, shared_ptr<WrapTFileOutput>> files_by_id;

    for(const auto& calibration : DataMap)
    {
        const string& calibrationID = calibration.first;
        if(calibrationID.empty()) {
            LOG(WARNING) << "Found empty calibrationID, ignored.";
            continue;
        }
        auto pos = calibrationID.find_first_of('/');
        string filename = folder + "/" +
                          calibrationID.substr(0, pos)
                          + ".root";
        auto it_file = files_by_filename.find(filename);
        if(it_file == files_by_filename.end()) {
            auto file = make_shared<WrapTFileOutput>(filename,
                                                     WrapTFileOutput::mode_t::recreate);

            files_by_filename.insert(make_pair(filename, file));
            files_by_id.insert(make_pair(calibrationID, file));
        }
        else {
            files_by_id.insert(make_pair(calibrationID, it_file->second));
        }

    }

    for(const auto& it_map : files_by_id)
    {
        WriteToTree(*it_map.second, *DataMap.find(it_map.first));
    }
}

uint32_t DataBase::getDepth(const TID& tid, const string& calibrationID) const
{
    uint32_t current_depth = 0;
    auto& calibPairs = DataMap.at(calibrationID);

    for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
    {
        interval<TID> idint(rit->FirstID,rit->LastID);
        if (idint.Contains(tid))
            return current_depth;
        current_depth++;
    }
    return current_depth;
}



uint32_t DataBase::GetNumberOfDataPoints(const string& calibrationID) const
{
    try
    {
        return DataMap.at(calibrationID).size();
    }
    catch (out_of_range)
    {
        return 0;
    }
}

const std::list<TID> DataBase::GetChangePoints(const string& calibrationID) const
{
    if ( DataMap.count(calibrationID) == 0)
        return {};

    uint32_t depth = 0;
    list<TID> ids;

    auto& calibPairs = DataMap.at(calibrationID);

    for(auto rit = calibPairs.rbegin(); rit != calibPairs.rend(); ++rit)
    {
        //changepoint is one after the last element;
        auto inclastID(rit->LastID);
        ++inclastID;

        if (isValid(rit->FirstID,calibrationID,depth) )
            ids.push_back(rit->FirstID);
        if (isValid(inclastID,calibrationID,depth) )
            ids.push_back(inclastID);

        depth++;
    }

    ids.sort();

    return ids;
}


