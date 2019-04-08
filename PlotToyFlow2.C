#include <TMath.h>
#include <TFile.h>

void PlotToyFlow2() {
    
    int i;
    double pi = TMath::Pi();
    const int n = 5;

    //TString sFileName = "toyFlow_1000000-events_random-psi.root";
    TString sFileName = "toyFlow_50000-events_dNeta-1500_random-psi.root";
    //TString sFileName = "toyFlow.root";
    //TString sFileName = "toyFlow_100000-events.root";
    TFile *fIn = TFile::Open(sFileName, "read");

    // phi, pT, eta, multiplicity
    TH1D *hPhi = (TH1D*)fIn->Get("hPhi");
    TH1D *hPt = (TH1D*)fIn->Get("hPt");
    TH1D *hEta = (TH1D*)fIn->Get("hEta");
    TH1D *hMultiplicity = (TH1D*)fIn->Get("hMultiplicity");

    // Resolution parameter
    TH1D *hRsub[n];
    for (i = 0; i < n; i++) {
        hRsub[i] = (TH1D*)fIn->Get(Form("hRsub%02i", i+1));
    }

    // For investigation of event plane resolution
    TH1D *hPsiA[n];
    TH1D *hPsiB[n];
    TH1D *hPsiDiff[n];
    TH1D *hPsiDiffTimesn[n];
    TH2D *hPsiAPsiB[n];
    for (i = 0; i < n; i++) {
        hPsiA[i] = (TH1D*)fIn->Get(Form("hPsiA%02i", i+1));
        hPsiB[i] = (TH1D*)fIn->Get(Form("hPsiB%02i", i+1));
        hPsiDiff[i] = (TH1D*)fIn->Get(Form("hPsiDiff%02i", i+1));
        hPsiDiffTimesn[i] = (TH1D*)fIn->Get(Form("hPsiDiffTimesn%02i", i+1));
        hPsiAPsiB[i] = (TH2D*)fIn->Get(Form("hPsiAPsiB%02i", i+1));
    }

    TH1D *hPsiANeg[n];
    TH1D *hPsiBNeg[n];
    TH1D *hPsiDiffNeg[n];
    TH1D *hPsiDiffTimesnNeg[n];
    TH2D *hPsiAPsiBNeg[n];
    for (i = 0; i < n; i++) {
        hPsiANeg[i] = (TH1D*)fIn->Get(Form("hPsiANeg%02i", i+1));
        hPsiBNeg[i] = (TH1D*)fIn->Get(Form("hPsiBNeg%02i", i+1));
        hPsiDiffNeg[i] = (TH1D*)fIn->Get(Form("hPsiDiffNeg%02i", i+1));
        hPsiDiffTimesnNeg[i] = (TH1D*)fIn->Get(Form("hPsiDiffTimesnNeg%02i", i+1));
        hPsiAPsiBNeg[i] = (TH2D*)fIn->Get(Form("hPsiAPsiBNeg%02i", i+1));
    }

    //=====PLOTS=====

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->Divide(2, 2);
    c1->cd(1); hPhi->Draw();
    c1->cd(2); hPt->Draw();
    c1->cd(3); hEta->Draw();
    c1->cd(4); hMultiplicity->Draw();

    TCanvas *c5 = new TCanvas("c5", "c5");
    c5->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c5->cd(i+1); hRsub[i]->Draw();
    }

    TCanvas *c7 = new TCanvas("c7", "c7");
    c7->Divide(3, 2);
    c7->cd(1); hPsiDiff[0]->Draw();
    c7->cd(2); hPsiDiffTimesn[0]->Draw();
    c7->cd(3); hPsiAPsiB[0]->Draw();
    c7->cd(4); hPsiDiffNeg[0]->Draw();
    c7->cd(5); hPsiDiffTimesnNeg[0]->Draw();
    c7->cd(6); hPsiAPsiBNeg[0]->Draw();

    TCanvas *c8 = new TCanvas("c8", "c8");
    c8->Divide(3, 2);
    c8->cd(1); hPsiDiff[1]->Draw();
    c8->cd(2); hPsiDiffTimesn[1]->Draw();
    c8->cd(3); hPsiAPsiB[1]->Draw();
    c8->cd(4); hPsiDiffNeg[1]->Draw();
    c8->cd(5); hPsiDiffTimesnNeg[1]->Draw();
    c8->cd(6); hPsiAPsiBNeg[1]->Draw();

    TCanvas *c9 = new TCanvas("c9", "c9");
    c9->Divide(3, 2);
    c9->cd(1); hPsiDiff[2]->Draw();
    c9->cd(2); hPsiDiffTimesn[2]->Draw();
    c9->cd(3); hPsiAPsiB[2]->Draw();
    c9->cd(4); hPsiDiffNeg[2]->Draw();
    c9->cd(5); hPsiDiffTimesnNeg[2]->Draw();
    c9->cd(6); hPsiAPsiBNeg[2]->Draw();

    TCanvas *c10 = new TCanvas("c10", "c10");
    c10->Divide(3, 2);
    c10->cd(1); hPsiDiff[3]->Draw();
    c10->cd(2); hPsiDiffTimesn[3]->Draw();
    c10->cd(3); hPsiAPsiB[3]->Draw();
    c10->cd(4); hPsiDiffNeg[3]->Draw();
    c10->cd(5); hPsiDiffTimesnNeg[3]->Draw();
    c10->cd(6); hPsiAPsiBNeg[3]->Draw();

    TCanvas *c11 = new TCanvas("c11", "c11");
    c11->Divide(3, 2);
    c11->cd(1); hPsiDiff[4]->Draw();
    c11->cd(2); hPsiDiffTimesn[4]->Draw();
    c11->cd(3); hPsiAPsiB[4]->Draw();
    c11->cd(4); hPsiDiffNeg[4]->Draw();
    c11->cd(5); hPsiDiffTimesnNeg[4]->Draw();
    c11->cd(6); hPsiAPsiBNeg[4]->Draw();
 }

