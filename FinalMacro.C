#include <TSystem.h>

void FinalMacro(bool bPlotCentData=true, TString sIn="toyFlow.root") {
    gROOT->ProcessLine(".I src/ResIter.h");

    if (gSystem->AccessPathName(sIn)) {
        cout << "File \"" << sIn << "\" not found!\n";
        return 0;
    }

    if (bPlotCentData) {
        cout << "Run MakeCentralityGraphs.C: " << Form(".x MakeCentralityGraphs.C(\"%s\")",sIn.Data()) << "\n";
        gROOT->ProcessLine(Form(".x MakeCentralityGraphs.C(\"%s\")",sIn.Data()));
        cout << "Run PlotCentralityData.C\n";
        gROOT->ProcessLine(".x PlotCentralityData.C");
    } else {
        cout << "Run MakeGraphs.C: " << Form(".x MakeGraphs.C(\"%s\")",sIn.Data()) << "\n";
        gROOT->ProcessLine(Form(".x MakeGraphs.C(\"%s\")",sIn.Data()));
        cout << "Run PlotVn.C\n";
        gROOT->ProcessLine(".x PlotVn.C");
    }

}
