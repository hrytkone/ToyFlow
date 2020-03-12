#define CENTDST_N 7

#include "src/JConst.h"

void checkUnderOverFlow( TH1 *h );
double GetRes(double rab, double rac, double rba);
double GetResError(double rab, double rabErr,
                   double rac, double racErr,
                   double rba, double rbaErr);
double GetVnError(double vobs, double vobsErr, double res, double resErr);

void PlotCentralityData3sub(TString sInputName = "input.root", TString sOutputName = "output.pdf", int n = 2, int iVnDet = 0) {

    TFile *fIn = TFile::Open(sInputName, "read");

    Float_t binLower[CENTDST_N+1] = {0, 5, 10, 20, 30, 40, 50, 60};

    double centvn[nCoef][CENTDST_N] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0277401, 0.04488324, 0.06521883, 0.08433443, 0.09597485, 0.10087206, 0.09925828},
        {0.02039728, 0.02369955, 0.02670301, 0.02950095, 0.03118808, 0.03120636, 0.02918556},
        {0.01013229, 0.01171893, 0.0131265, 0.01479335, 0.0159713, 0.01644628, 0.01535014},
        {0.00415816, 0.00467961, 0.00528238, 0.006501, 0.0068885, 0.00690379, 0.00575251}
    };

    TH1D *hInput = new TH1D("hInput", "hInput", CENTDST_N, binLower);
    for (int i=1; i<=CENTDST_N; i++) {
        hInput->SetBinContent(i, centvn[n-1][i-1]);
    }

    TH1D *hVnObs[CENTBINS_N];
    TH1D *hRsubAB[CENTBINS_N];
    TH1D *hRsubAC[CENTBINS_N];
    TH1D *hRsubBA[CENTBINS_N];
    for (int i=0; i<CENTBINS_N; i++) {
        hVnObs[i] = (TH1D*)fIn->Get(Form("hVnObsH%02iD%02iCENT%02i", n, iVnDet, i));
        checkUnderOverFlow(hVnObs[i]);
        hRsubAB[i] = (TH1D*)fIn->Get(Form("hRsubAB%dCENT%02d", n, i));
        checkUnderOverFlow(hRsubAB[i]);
        hRsubAC[i] = (TH1D*)fIn->Get(Form("hRsubAC%dCENT%02d", n, i));
        checkUnderOverFlow(hRsubAC[i]);
        hRsubBA[i] = (TH1D*)fIn->Get(Form("hRsubBA%dCENT%02d", n, i));
        checkUnderOverFlow(hRsubBA[i]);
    }

    double r[CENTBINS_N], rerr[CENTBINS_N];
    double vobs[CENTBINS_N], vobserr[CENTBINS_N];
    double v[CENTBINS_N], verr[CENTBINS_N];
    for (int i=0; i<CENTBINS_N; i++) {
        cout << "\nCENTRALITY BIN " << i << " : \n";

        double res = GetRes(hRsubAB[i]->GetMean(), hRsubAC[i]->GetMean(), hRsubBA[i]->GetMean());
        if (isnan(res)) res = 0.0;

        r[i] = res;
        vobs[i] = hVnObs[i]->GetMean();
        if (res==0.0) {
            v[i] = 0.0;
        } else {
            v[i] = vobs[i]/res;
        }

        rerr[i] = GetResError(hRsubAB[i]->GetMean(), hRsubAB[i]->GetMeanError(),
                              hRsubAC[i]->GetMean(), hRsubAC[i]->GetMeanError(),
                              hRsubBA[i]->GetMean(), hRsubBA[i]->GetMeanError());
        vobserr[i] = hVnObs[i]->GetMeanError();
        verr[i] = GetVnError(vobs[i], vobserr[i], r[i], rerr[i]);

        cout << "\n   R" << nCoef << " : " << res << " +- " << rerr[i] << endl;
        cout << "   vobs" << nCoef << " : " << vobs[i] << " +- " << vobserr[i] << endl;
        cout << "   v" << nCoef << " : " << v[i] << " +- " << verr[i] << endl;
    }

    double mrkSize = 1.0;
    int mrkStyle = 20;

    TGraphErrors *gR = new TGraphErrors(CENTDST_N);
    TGraphErrors *gVobs = new TGraphErrors(CENTDST_N);
    TGraphErrors *gV = new TGraphErrors(CENTDST_N);

    for(int i=0; i<CENTDST_N; i++) {
        if (i==0 || i==1) {
            gR->SetPoint(i, double(i+1)*5.0-2.5, r[i]);
            gR->SetPointError(i, 0.0, rerr[i]);

            gVobs->SetPoint(i, double(i+1)*5.0-2.5, vobs[i]);
            gVobs->SetPointError(i, 0.0, vobserr[i]);

            gV->SetPoint(i, double(i+1)*5.0-2.5, v[i]);
            gV->SetPointError(i, 0.0, verr[i]);
        } else {
            gR->SetPoint(i, double(i)*10.0-5.0, r[i]);
            gR->SetPointError(i, 0.0, rerr[i]);

            gVobs->SetPoint(i, double(i)*10.0-5.0, vobs[i]);
            gVobs->SetPointError(i, 0.0, vobserr[i]);

            gV->SetPoint(i, double(i)*10.0-5.0, v[i]);
            gV->SetPointError(i, 0.0, verr[i]);
        }
    }

    TLegend *leg = new TLegend(0.15, 0.7, 0.3, 0.85);
    leg->AddEntry(hInput , "Input", "l");
    leg->AddEntry(gV , Form("v_{%01i}", 2), "p");
    leg->AddEntry(gVobs , Form("v_{obs%01i}", 2), "p");

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 500);
    //TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
    c1->Divide(2, 1);

    c1->cd(1);
    gR->Draw("AP");
    gR->GetXaxis()->CenterTitle(true);
    gR->GetYaxis()->CenterTitle(true);
    gR->GetYaxis()->SetRangeUser(0.5,1.0);
    gR->SetTitle(Form("Resolution for V0C; centrality; R_{%01i}", 2));
    gR->SetMarkerColor(4);
    gR->SetMarkerStyle(mrkStyle);
    gR->SetMarkerSize(mrkSize);
    gR->Draw("AP");
    c1->Update();

    c1->cd(2);
    gV->Draw("AP");
    gV->GetXaxis()->CenterTitle(true);
    gV->GetYaxis()->CenterTitle(true);
    gV->GetYaxis()->SetRangeUser(0.0, 0.25);
    gV->SetTitle(Form("v_{%01i} for TPC; centrality; v_{%01i}", 2, 2));
    gV->SetMarkerColor(4);
    gV->SetMarkerStyle(mrkStyle);
    gV->SetMarkerSize(mrkSize);
    gV->Draw("AP");
    c1->Update();

    gVobs->Draw("P SAME");
    gVobs->SetMarkerColor(2);
    gVobs->GetYaxis()->SetRangeUser(0.0, 0.45);
    gVobs->SetMarkerStyle(mrkStyle);
    gVobs->SetMarkerSize(mrkSize);

    hInput->Draw("HIST SAME");
    leg->Draw();

    c1->SaveAs(sOutputName);
}

//______________________________________________________________________________

void checkUnderOverFlow( TH1 *h ) {
        if(h->GetBinContent(0)>0) cout << h->GetName() << " underflow bin not empty: " << h->GetBinContent(0) << endl;
        if(h->GetBinContent(h->GetXaxis()->GetNbins()+1)>0) cout << h->GetName() << " overflow bin not empty: " << h->GetBinContent(h->GetXaxis()->GetNbins()+1) << endl;
}

double GetRes(double rab, double rac, double rba) {
    return TMath::Sqrt((rab * rac)/rba);
}

double GetResError(double rab, double rabErr,
                   double rac, double racErr,
                   double rba, double rbaErr) {
    return 0.5*GetRes(rab, rac, rba)*TMath::Sqrt(TMath::Power(rabErr/rab, 2.0) +
    TMath::Power(racErr/rac, 2.0) + TMath::Power(rbaErr/rba, 2.0));
}

double GetVnError(double vobs, double vobsErr, double res, double resErr) {
    return (vobs/res)*TMath::Sqrt(TMath::Power(vobsErr/vobs, 2.0) +
    TMath::Power(resErr/res, 2.0));
}
