#include <TMath.h>
#include <TFile.h>

double ShiftCorrection(double x, double correction);
double TwistCorrectionX(double x, double y, double lambdaMinus, double lambdaPlus);
double TwistCorrectionY(double y, double x, double lambdaMinus, double lambdaPlus);
double RescalingCorrectionX(double x, double aPlus);
double RescalingCorrectionY(double y, double aMinus);
void Shift(TH1D *hX, TH1D *hY, TH1D *hXshifted, TH1D *hYshifted, double cm, double sm);
void Twist(TH1D *hX, TH1D *hY, TH1D *hXtwisted, TH1D *hYtwisted, double lambdaMinus, double lambdaPlus);
void Rescaling(TH1D *hX, TH1D *hY, TH1D *hXrescaled, TH1D *hYrescaled, double aMinus, double aPlus);

void PlotToyFlowCorrections() {

    int i;
    const int n = 5;

    TString sFileName = "toyFlow_10000-events_60_pi4.root";
    TFile *fIn = TFile::Open(sFileName, "read");

    TH1D *hPhi = (TH1D*)fIn->Get("hPhi");
    TH1D *hPhiNonuni = (TH1D*)fIn->Get("hPhiNonuni");

    // Q-vector for uniform phi-distribution
    TH1D *hQx[n];
    TH1D *hQy[n];
    for (i=0; i<n; i++) {
        hQx[i] = (TH1D*)fIn->Get(Form("hQx%02i", i+1));
        hQy[i] = (TH1D*)fIn->Get(Form("hQy%02i", i+1));
    }

    // Q-vector and cos(m*phi), sin(m*phi) for non-uniform phi-distribution
    TH1D *hQxNonuni[n];
    TH1D *hQyNonuni[n];
    TH1D *hCosPhi[n];
    TH1D *hSinPhi[n];
    TH1D *hCosPhi2[n];
    TH1D *hSinPhi2[n];
    for (i=0; i<n; i++) {
        hQxNonuni[i] = (TH1D*)fIn->Get(Form("hQxNonuni%02i", i+1));
        hQyNonuni[i] = (TH1D*)fIn->Get(Form("hQyNonuni%02i", i+1));
        hCosPhi[i] = (TH1D*)fIn->Get(Form("hCosPhi%02i", i+1));
        hSinPhi[i] = (TH1D*)fIn->Get(Form("hSinPhi%02i", i+1));
        hCosPhi2[i] = (TH1D*)fIn->Get(Form("hCosPhi2%02i", i+1));
        hSinPhi2[i] = (TH1D*)fIn->Get(Form("hSinPhi2%02i", i+1));
    }

    int k = 1;

    double cn = hCosPhi[k]->GetMean();
    double sn = hSinPhi[k]->GetMean();

    cout << "cn=" << cn << endl;
    cout << "sn=" << sn << endl;
    cout << "" << endl;

    double cn2 = hCosPhi2[k]->GetMean();
    double sn2 = hSinPhi2[k]->GetMean();

    cout << "cn2=" << cn2 << endl;
    cout << "sn2=" << sn2 << endl;
    cout << "" << endl;

    double aMinus = 1 - cn2;
    double aPlus = 1 + cn2;

    cout << "aMinus=" << aMinus << endl;
    cout << "aPlus=" << aPlus << endl;
    cout << "" << endl;

    double lambdaMinus = sn2/aMinus;
    double lambdaPlus = sn2/aPlus;

    cout << "lambdaMinus=" << lambdaMinus << endl;
    cout << "lambdaPlus=" << lambdaPlus << endl;
    cout << "" << endl;

    // Corrections to mean values
    double QxShifted = ShiftCorrection(hQxNonuni[k]->GetMean(), cn);
    double QyShifted = ShiftCorrection(hQyNonuni[k]->GetMean(), sn);

    double QxTwisted = TwistCorrectionX(QxShifted, QyShifted, lambdaMinus, lambdaPlus);
    double QyTwisted = TwistCorrectionY(QxShifted, QyShifted, lambdaMinus, lambdaPlus);

    double QxRescaled = RescalingCorrectionX(QxTwisted, aPlus);
    double QyRescaled = RescalingCorrectionY(QyTwisted, aMinus);

    cout << "Uniform" << endl;
    cout << "===========" << endl;
    cout << "  Qx=" << hQx[k]->GetMean() << endl;
    cout << "  Qy=" << hQy[k]->GetMean() << endl;
    cout << "" << endl;
    cout << "Non-uniform" << endl;
    cout << "===========" << endl;
    cout << "  Non-corrected:" << endl;
    cout << "    Qx=" << hQxNonuni[k]->GetMean() << endl;
    cout << "    Qy=" << hQyNonuni[k]->GetMean() << endl;
    cout << "  Shifted:" << endl;
    cout << "    Qx=" << QxShifted << endl;
    cout << "    Qy=" << QyShifted << endl;
    cout << "  Twisted:" << endl;
    cout << "    Qx=" << QxTwisted << endl;
    cout << "    Qy=" << QyTwisted << endl;
    cout << "  Rescaled:" << endl;
    cout << "    Qx=" << QxRescaled << endl;
    cout << "    Qy=" << QyRescaled << endl;

    // Corrections to histograms
    TH1D *hXshifted = new TH1D("hXshifted", "Qx shifted", 401, -25.0, 25.0);
    TH1D *hYshifted = new TH1D("hYshifted", "Qy shifted", 401, -25.0, 25.0);
    TH1D *hXtwisted = new TH1D("hXtwisted", "Qx twisted", 401, -25.0, 25.0);
    TH1D *hYtwisted = new TH1D("hYtwisted", "Qy twisted", 401, -25.0, 25.0);
    TH1D *hXrescaled = new TH1D("hXrescaled", "Qx rescaled", 401, -25.0, 25.0);
    TH1D *hYrescaled = new TH1D("hYrescaled", "Qy rescaled", 401, -25.0, 25.0);

    Shift(hQxNonuni[k], hQyNonuni[k], hXshifted, hYshifted, cn, sn);
    Twist(hXshifted, hYshifted, hXtwisted, hYtwisted, lambdaMinus, lambdaPlus);
    Rescaling(hXtwisted, hYtwisted, hXrescaled, hYrescaled, aMinus, aPlus);

    // Plot everything
    hQx[k]->SetLineColor(kRed);
    hQxNonuni[k]->SetLineColor(kGreen);
    hXshifted->SetLineColor(kCyan);
    hXtwisted->SetLineColor(kMagenta);
    hXrescaled->SetLineColor(kBlue);

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->cd();
    hQx[k]->Draw("hist");
    hQxNonuni[k]->Draw("hist same");
    hXshifted->Draw("hist same");
    hXtwisted->Draw("hist same");
    hXrescaled->Draw("hist same");

    TLegend *legend1 = new TLegend(0.1,0.7,0.48,0.9);
    legend1->SetHeader("Corrections for Qx", "C"); // option "C" allows to center the header
    legend1->AddEntry(hQx[k], "Uniform", "l");
    legend1->AddEntry(hQxNonuni[k], "Nonuniform", "l");
    legend1->AddEntry(hXshifted, "Shift", "l");
    legend1->AddEntry(hXtwisted, "Twist", "l");
    legend1->AddEntry(hXrescaled, "Rescaling", "l");
    legend1->Draw();

    hQy[k]->SetLineColor(kRed);
    hQyNonuni[k]->SetLineColor(kGreen);
    hYshifted->SetLineColor(kCyan);
    hYtwisted->SetLineColor(kMagenta);
    hYrescaled->SetLineColor(kBlue);

    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->cd();
    hQy[k]->Draw("hist");
    hQyNonuni[k]->Draw("hist same");
    hYshifted->Draw("hist same");
    hYtwisted->Draw("hist same");
    hYrescaled->Draw("hist same");

    TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
    legend2->SetHeader("Corrections for Qy", "C"); // option "C" allows to center the header
    legend2->AddEntry(hQy[k], "Uniform", "l");
    legend2->AddEntry(hQyNonuni[k], "Nonuniform", "l");
    legend2->AddEntry(hYshifted, "Shift", "l");
    legend2->AddEntry(hYtwisted, "Twist", "l");
    legend2->AddEntry(hYrescaled, "Rescaling", "l");
    legend2->Draw();

}

double ShiftCorrection(double x, double correction) {
    return x - correction;
}

double TwistCorrectionX(double x, double y, double lambdaMinus, double lambdaPlus) {
    return (x - lambdaMinus*y)/(1-lambdaMinus*lambdaPlus);
}

double TwistCorrectionY(double y, double x, double lambdaMinus, double lambdaPlus) {
    return (y - lambdaPlus*x)/(1-lambdaMinus*lambdaPlus);
}

double RescalingCorrectionX(double x, double aPlus) {
    return x/aPlus;
}

double RescalingCorrectionY(double y, double aMinus) {
    return y/aMinus;
}

void Shift(TH1D *hX, TH1D *hY, TH1D *hXshifted, TH1D *hYshifted, double cm, double sm) {
    int i, j;
    double xCont, yCont, binCenter;

    for (i=0; i<hX->GetNbinsX(); i++) {
        xCont = hX->GetBinContent(i);
        yCont = hY->GetBinContent(i);

        if (xCont>0) {
            for (j=0; j<xCont; j++) {
                binCenter = hXshifted->GetXaxis()->GetBinCenter(i);
                hXshifted->Fill(binCenter-cm);
            }
        }

        if (yCont>0) {
            for (j=0; j<yCont; j++) {
                binCenter = hYshifted->GetXaxis()->GetBinCenter(i);
                hYshifted->Fill(binCenter-sm);
            }
        }
    }
}

void Twist(TH1D *hX, TH1D *hY, TH1D *hXtwisted, TH1D *hYtwisted, double lambdaMinus, double lambdaPlus) {
    int i, j;
    double xCont, yCont, xBinCenter, yBinCenter;

    for (i=0; i<hX->GetNbinsX(); i++) {
        xCont = hX->GetBinContent(i);
        yCont = hY->GetBinContent(i);

        if (xCont>0) {
            for (j=0; j<xCont; j++) {
                xBinCenter = hXtwisted->GetXaxis()->GetBinCenter(i);
                yBinCenter = hYtwisted->GetXaxis()->GetBinCenter(i);
                hXtwisted->Fill((xBinCenter-lambdaMinus*yBinCenter)/(1-lambdaMinus*lambdaPlus));
            }
        }

        if (xCont>0) {
            for (j=0; j<yCont; j++) {
                xBinCenter = hXtwisted->GetXaxis()->GetBinCenter(i);
                yBinCenter = hYtwisted->GetXaxis()->GetBinCenter(i);
                hYtwisted->Fill((yBinCenter-lambdaPlus*xBinCenter)/(1-lambdaMinus*lambdaPlus));
            }
        }
    }
}

void Rescaling(TH1D *hX, TH1D *hY, TH1D *hXrescaled, TH1D *hYrescaled, double aMinus, double aPlus) {
    int i, j;
    double xCont, yCont, xBinCenter, yBinCenter;

    for (i=0; i<hX->GetNbinsX(); i++) {
        xCont = hX->GetBinContent(i);
        yCont = hY->GetBinContent(i);

        if (xCont>0) {
            for (j=0; j<xCont; j++) {
                xBinCenter = hXrescaled->GetXaxis()->GetBinCenter(i);
                hXrescaled->Fill(xBinCenter/aPlus);
            }
        }

        if (yCont>0) {
            for (j=0; j<xCont; j++) {
                yBinCenter = hYrescaled->GetXaxis()->GetBinCenter(i);
                hYrescaled->Fill(yBinCenter/aMinus);
            }
        }
    }
}
