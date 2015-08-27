#pragma once

#include "Manager_traits.h"

#include "TGFrame.h"
#include "KeySymbols.h"

#include <functional>
#include <list>
#include <map>

class TGStatusBar;
class TGWidget;

namespace ant {
namespace calibration {
namespace gui {

class CalCanvas;
class Manager;
class ProgressBar;

struct ManagerWindowMode {
    ManagerWindowMode() :
        gotoNextSlice(true),
        autoContinue(true),
        showEachFit(true),
        channelStep(1),
        requestChannel(-1)
    {}

    bool gotoNextSlice;
    bool autoContinue;
    bool showEachFit;
    int  channelStep;
    int  requestChannel;
};

class ManagerWindow : public TGMainFrame, public ManagerWindow_traits
{
private:
    std::list<CalCanvas*> canvases;
    TGHorizontalFrame* frame_canvases = nullptr;

    TGStatusBar* statusbar = nullptr;
    ProgressBar* progress_channel = nullptr;
    ProgressBar* progress_slice = nullptr;

    std::list<TGWidget*> nonfinish_widgets;

    std::map<EKeySym, TGTextButton*> keys;
    Manager* manager = nullptr;


    void CreateToolbar(TGVerticalFrame* frame);
    void UpdateLayout();
    void RunManager();

public:
    ManagerWindow(Manager* manager_);
    virtual Bool_t HandleKey(Event_t *event) override;
    virtual ~ManagerWindow();
    virtual gui::CalCanvas* AddCalCanvas(const std::string& name = "") override;

    ManagerWindowMode Mode;

    void SetProgressMax(unsigned slices, unsigned channels);
    void SetProgress(unsigned slice, unsigned channel);

    void SetFinishMode(bool flag);

    ManagerWindow(const ManagerWindow&) = delete;
    ManagerWindow& operator=(const ManagerWindow&) = delete;
};
}
}
}
