void FinalMacro(TString sIn="toyFlow.root") {
    gROOT->ProcessLine(".I /home/heimarry/alice/sw/slc7_x86-64/boost/v1.59.0-1/include/");
    gROOT->ProcessLine(".I src/ResIter.h");

    cout << "Run MakeGraphs.C: " << Form(".x MakeGraphs.C(\"%s\")",sIn.Data()) << endl;
    gROOT->ProcessLine(Form(".x MakeGraphs.C(\"%s\")",sIn.Data()));
    cout << "Run PlotVn.C" << endl;
    gROOT->ProcessLine(".x PlotVn.C");
    //gROOT->ProcessLine(".x TestIter.C");
}
