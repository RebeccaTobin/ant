#pragma once

//ant
#include "tree/TCalibrationData.h"
#include "tree/TDataRecord.h"
#include "DataBase.h"

//std
#include <list>
#include <string>
#include <memory>


namespace ant
{
namespace calibration
{

class DataAccess
{
public:
   /**
    *  \brief GetData Query the calibration database for specific TID
    *  \param calibrationID Calibration ID
    *  \param eventID event ID
    *  \param cdata   Reference to a TCalibrationData, data will be writter here
    *  \return true if valid data was found
    */
    virtual void Add(const TCalibrationData& data) = 0;

    virtual bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) = 0;

    virtual const std::list<TID> GetChangePoints(const std::string& calibrationID) = 0;

    virtual bool GetLastEntry(const std::string& calibrationID, TCalibrationData& cdata) = 0;
};


class DataManager: public DataAccess
{

private:

    std::string calibrationDataFolder;
    bool changedDataBase;
    std::unique_ptr<DataBase> dataBase;

    void Init();



public:
    DataManager(const std::string& calibrationDataFolder_):
        calibrationDataFolder(calibrationDataFolder_),
        changedDataBase(false)
    {}


    ~DataManager()
    {
        if (changedDataBase)
            dataBase->WriteToFolder(calibrationDataFolder);
    }


    void Add(const TCalibrationData& data) override;

    bool GetData(const std::string& calibrationID, const TID& eventID, TCalibrationData& cdata) override;

    const std::list<TID> GetChangePoints(const std::string& calibrationID) override;

    std::uint32_t GetNumberOfCalibrations();

    std::uint32_t GetNumberOfDataPoints(const std::string& calibrationID);


    bool GetLastEntry(const std::string& calibrationID, TCalibrationData& cdata) override;

};


} //namespace calibration
} //namespace ant
