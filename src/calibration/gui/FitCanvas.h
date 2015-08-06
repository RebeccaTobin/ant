#pragma once

#include "calibration/gui/GUIbase.h"
#include "calibration/gui/FitFunction.h"

#include <list>
#include <memory>
#include <string>
#include <stack>

#include "TCanvas.h"


class TH1;
class TRootCanvas;

namespace ant {
namespace calibration {
namespace gui {



class CalCanvas : public TCanvas, public update_notify_traits {
protected:
    Viewport getViewport();
    void ClearInidators();
    void SetupGUI();

    GUIIndicator* MakeVerticalIndicatorLine(VirtualKnob& knob);
    GUIIndicator* MakeHorizontalIndicatorLine(VirtualKnob& knob);
    GUIIndicator* MakeGUIElement(VirtualKnob& knob);

    std::list<GUIIndicator*> indicators;
    FitFunction* func = nullptr;
    TH1* hist = nullptr;

    virtual void HandleKeypress(const char key);

    std::stack<FitFunction::SavedState_t> UndoStack;

    TRootCanvas* rootcanvas;

public:
    CalCanvas(const std::string& name);
    virtual ~CalCanvas();

    virtual void Show(TH1* h, FitFunction* f);

    virtual void update_me() override;

    virtual void ShowGuidelines(TObject*, const Int_t, const char, const bool) override;

    virtual void Fit();

    virtual void UndoPush();
    virtual void UndoPop();

    virtual void ConnectReturnFunc(const char* receiver_class, void* receiver, const char* slot);

    virtual void Execute(const char *method, const char *params, Int_t *error);

    virtual void Update() override;




    /**
     * @brief HandleInput: Override default to catch keyboard inputs
     * @param button
     * @param x
     * @param y
     */
    virtual void HandleInput(EEventType button, Int_t x, Int_t y) override;

};
}
}
}
