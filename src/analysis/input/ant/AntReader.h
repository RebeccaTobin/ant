#pragma once

#include "analysis/input/DataReader.h"

#include "unpacker/Unpacker.h"
#include "tree/UnpackerWriter.h"

#include "reconstruct/Reconstruct_traits.h"

#include "Rtypes.h"

#include <memory>
#include <string>


class TTree;

namespace ant {

class TEvent;

namespace analysis {

namespace data {
    class Event;
}

namespace input {


class AntReader : public DataReader {
protected:
    std::unique_ptr<Unpacker::Reader> reader;
    std::unique_ptr<tree::UnpackerWriter> writer;
    std::unique_ptr<Reconstruct_traits> reconstruct;
    bool haveReconstruct;


    bool writeUncalibrated;
    bool writeCalibrated;

    std::unique_ptr<TSlowControl> buffered_slowcontrol;

public:
    AntReader(std::unique_ptr<Unpacker::Reader> unpacker_reader,
                      std::unique_ptr<Reconstruct_traits> reconstruct = nullptr);
    virtual ~AntReader();
    AntReader(const AntReader&) = delete;
    AntReader& operator= (const AntReader&) = delete;

    void EnableUnpackerWriter(const std::string& outputfile,
                              bool uncalibratedDetectorReads = false,
                              bool calibratedDetectorReads = false);

    // DataReader interface
    virtual bool IsSource() override { return true; }
    virtual bool ReadNextEvent(data::Event& event) override;
    virtual std::unique_ptr<TSlowControl> ReadNextSlowControl() override;

    double PercentDone() const override;
};

}
}
}
