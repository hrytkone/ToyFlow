#include <TMath.h>
#include <TFile.h>

#include "ResIter.h"

void PlotToyFlow() {
    
    int i;
    double pi = TMath::Pi();
    const int n = 5;

    //TString sFileName = "toyFlow_50000-events_dNeta-1000_random-psi.root";
    TString sFileName = "toyFlow.root";
    TFile *fIn = TFile::Open(sFileName, "read");

    // phi, pT, eta, multiplicity
    TH1D *hPhi = (TH1D*)fIn->Get("hPhi");
    TH1D *hPt = (TH1D*)fIn->Get("hPt");
    TH1D *hEta = (TH1D*)fIn->Get("hEta");
    TH1D *hMultiplicity = (TH1D*)fIn->Get("hMultiplicity");

    // Q-vector
    /**TH1D *hQx[n];
    TH1D *hQy[n];
    for (i = 0; i < n; i++) {
        hQx[i] = (TH1D*)fIn->Get(Form("hQx%02i", i+1));
        hQy[i] = (TH1D*)fIn->Get(Form("hQy%02i", i+1));
    }

    // observed vn
    TH1D *hVnObs[n];
    for (i = 0; i < n; i++) {
        hVnObs[i] = (TH1D*)fIn->Get(Form("hVnObs%02i", i+1));
    }**/

    // Resolution parameter
    TH1D *hRsub[n];
    for (i = 0; i < n; i++) {
        hRsub[i] = (TH1D*)fIn->Get(Form("hRsub%02i", i+1));
    }

    // SP methods
    TH1D *hSPnom[n];
    TH1D *hSPdenom[n];
    for (i = 0; i < n; i++) {
        hSPnom[i] = (TH1D*)fIn->Get(Form("hSPnomV%02i", i+1));
        hSPdenom[i] = (TH1D*)fIn->Get(Form("hSPdenomV%02i", i+1));
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

    //=====vn{EP}=====

    /**double vnEP[n];
    double Rinit;

    double ksi0 = 2.0;
    double err = 0.0001;
    double ksi[n];
    double R[n];

    for (i = 0; i < n; i++) {
        vnEP[i] = hVnObs[i]->GetMean();
        Rinit = hRsub[i]->GetMean();
        ksi[i] = RkIter(ksi0, Rinit, n, err);
        R[i] = Rk(TMath::Sqrt(2)*ksi[i], n);
        vnEP[i] /= R[i];
    }
    
    cout << "   ===Rk===" << endl;
    for (i = 0; i < n; i++)
        cout << "   Rk" << i+1 << "=" << R[i] << endl;

    cout << "" << endl;
    cout << "   ===vn{EP}===" << endl;
    for (i = 0; i < n; i++)
        cout << "   v" << i+1 << "=" << vnEP[i] << endl;**/

    //=====vn{SP}=====
    
    double vnSP[n];

    for (i = 1; i < n; i++)
        vnSP[i] = hSPnom[i]->GetMean()/(2*TMath::Sqrt(200*hSPdenom[i]->GetMean()));

    cout << "" << endl;
    cout << "   ===vn{SP}===" << endl;
    for (i = 0; i < n; i++)
        cout << "   v" << i+1 << "=" << vnSP[i] << endl;

    //=====PLOTS=====

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->Divide(2, 2);
    c1->cd(1); hPhi->Draw();
    c1->cd(2); hPt->Draw();
    c1->cd(3); hEta->Draw();
    c1->cd(4); hMultiplicity->Draw();

    /**TCanvas *c2 = new TCanvas("c2", "c2");
    c2->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c2->cd(i+1); hQx[i]->Draw();
    }

    TCanvas *c3 = new TCanvas("c3", "c3");
    c3->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c3->cd(i+1); hQy[i]->Draw();
    }

    TCanvas *c4 = new TCanvas("c4", "c4");
    c4->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c4->cd(i+1); hVnObs[i]->Draw();
    }**/

    TCanvas *c5 = new TCanvas("c5", "c5");
    c5->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c5->cd(i+1); hRsub[i]->Draw();
    }

    /**gStyle->SetOptTitle(kFALSE);

    double x[5] = {1, 2, 3, 4, 5};
    TGraph *gr1 = new TGraph(n, x, vnEP);
    TGraph *gr2 = new TGraph(n, x, vnSP);

    TCanvas *c6 = new TCanvas("c6", "c6");

    gr1->SetTitle("vn{EP}");
    gr2->SetTitle("vn{SP}");

    gr1->SetLineWidth(0); 
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColor(46);
    gr1->Draw();

    gr2->SetLineWidth(0); 
    gr2->SetMarkerStyle(22);
    gr2->SetMarkerColor(9);
    gr2->Draw("P");

    gPad->BuildLegend();**/

    TCanvas *c7 = new TCanvas("c7", "c7");
    c7->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c7->cd(i+1); hPsiA[i]->Draw();
    }

    TCanvas *c8 = new TCanvas("c8", "c8");
    c8->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c8->cd(i+1); hPsiB[i]->Draw();
    }

    TCanvas *c9 = new TCanvas("c9", "c9");
    c9->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c9->cd(i+1); hPsiDiff[i]->Draw();
    }

    TCanvas *c10 = new TCanvas("c10", "c10");
    c10->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c10->cd(i+1); hPsiDiffTimesn[i]->Draw();
    }

    TCanvas *c11 = new TCanvas("c11", "c11");
    c11->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c11->cd(i+1); hPsiAPsiB[i]->Draw();
    }

    TCanvas *c12 = new TCanvas("c12", "c12");
    c12->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c12->cd(i+1); hPsiANeg[i]->Draw();
    }

    TCanvas *c13 = new TCanvas("c13", "c13");
    c13->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c13->cd(i+1); hPsiBNeg[i]->Draw();
    }

    TCanvas *c14 = new TCanvas("c14", "c14");
    c14->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c14->cd(i+1); hPsiDiffNeg[i]->Draw();
    }

    TCanvas *c15 = new TCanvas("c15", "c15");
    c15->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c15->cd(i+1); hPsiDiffTimesnNeg[i]->Draw();
    }

    TCanvas *c16 = new TCanvas("c16", "c16");
    c16->Divide(3, 2);
    for (i = 0; i < n; i++) {
        c16->cd(i+1); hPsiAPsiBNeg[i]->Draw();
    }
}

