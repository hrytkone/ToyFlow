#define CENTDST_N 7

#include "src/JConst.h"

void checkUnderOverFlow( TH1 *h );
double GetRes(double rab, double rac, double rba);
double GetResError(double rab, double rabErr,
                   double rac, double racErr,
                   double rba, double rbaErr);
double GetVnError(double vobs, double vobsErr, double res, double resErr);

const int nFiles = 3;
const int nFilesRead = 3;
const int nRef = 0;
TString fileName[nFiles] = {
    //"/home/heimarry/Desktop/toyflow-vntests/outputs-gran-on/toyFlow_20200312_PtDep0_Gran1_Scale1.0/data/toyFlow_20200312_PtDep0_Gran1_Scale1.0.root",
    //"/home/heimarry/Desktop/toyflow-vntests/outputs-gran-on/toyFlow_20200312_PtDep0_Gran1_Scale0.8/data/toyFlow_20200312_PtDep0_Gran1_Scale0.8.root",
    //"/home/heimarry/Desktop/toyflow-vntests/outputs-gran-on/toyFlow_20200312_PtDep0_Gran1_Scale0.65/data/toyFlow_20200312_PtDep0_Gran1_Scale0.65.root"
    //"output/toyFlow_20200429_PtDep0_Gran0_Scale1.00/toyFlow_20200429_PtDep0_Gran0_Scale1.00.root",
    //"output/toyFlow_20200429_PtDep0_Gran0_Scale0.80/toyFlow_20200429_PtDep0_Gran0_Scale0.80.root",
    //"output/toyFlow_20200429_PtDep0_Gran0_Scale0.65/toyFlow_20200429_PtDep0_Gran0_Scale0.65.root"
    "output/toyFlow_20200430_70perMulti_PtDep0_Gran0_Scale1.00/toyFlow_20200430_70perMulti_PtDep0_Gran0_Scale1.00.root",
    "output/toyFlow_20200430_70perMulti_PtDep0_Gran0_Scale0.80/toyFlow_20200430_70perMulti_PtDep0_Gran0_Scale0.80.root",
    "output/toyFlow_20200430_70perMulti_PtDep0_Gran0_Scale0.65/toyFlow_20200430_70perMulti_PtDep0_Gran0_Scale0.65.root",
};
TString sSame[nFiles] = {"AP", "SAME P", "SAME P"};
int gColor[nFiles] = {1,2,4};

int mMarker = 20;
double mSize = 1.0;

TGraphErrors *gR[nFiles];
TGraphErrors *gVobs[nFiles];
TGraphErrors *gV[nFiles];

TFile *fIn[nFiles];

void PlotCentralityData3sub(TString sOutputName = "output.pdf", int n = 2, int iVnDet = 0) {

    for(int iFil=0; iFil<nFiles; iFil++)
        fIn[iFil] = TFile::Open(fileName[iFil], "read");

    Float_t binLower[CENTDST_N+1] = {0, 5, 10, 20, 30, 40, 50, 60};

    double centvn[nCoef][CENTDST_N] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0277401, 0.04488324, 0.06521883, 0.08433443, 0.09597485, 0.10087206, 0.09925828},
        {0.02039728, 0.02369955, 0.02670301, 0.02950095, 0.03118808, 0.03120636, 0.02918556},
        {0.01013229, 0.01171893, 0.0131265, 0.01479335, 0.0159713, 0.01644628, 0.01535014},
        {0.00415816, 0.00467961, 0.00528238, 0.006501, 0.0068885, 0.00690379, 0.00575251}
    };

    double scale = 1.0;
    TH1D *hInput = new TH1D("hInput", "hInput", CENTDST_N, binLower);
    for (int i=1; i<=CENTDST_N; i++) {
        hInput->SetBinContent(i, scale * centvn[n-1][i-1]);
    }

    TH1D *hVnObs[CENTBINS_N];
    TH1D *hRsubAB[CENTBINS_N];
    TH1D *hRsubAC[CENTBINS_N];
    TH1D *hRsubBC[CENTBINS_N];

    for(int iFil=0; iFil<nFilesRead; iFil++) {

        gR[iFil] = new TGraphErrors(CENTDST_N);
        gVobs[iFil] = new TGraphErrors(CENTDST_N);
        gV[iFil] = new TGraphErrors(CENTDST_N);

        for (int i=0; i<CENTBINS_N; i++) {
            hVnObs[i] = (TH1D*)fIn[iFil]->Get(Form("hVnObsH%02iD%02iCENT%02i", n, iVnDet, i));
            checkUnderOverFlow(hVnObs[i]);
            hRsubAB[i] = (TH1D*)fIn[iFil]->Get(Form("hRsubAB%dCENT%02d", n, i));
            checkUnderOverFlow(hRsubAB[i]);
            hRsubAC[i] = (TH1D*)fIn[iFil]->Get(Form("hRsubAC%dCENT%02d", n, i));
            checkUnderOverFlow(hRsubAC[i]);
            hRsubBC[i] = (TH1D*)fIn[iFil]->Get(Form("hRsubBC%dCENT%02d", n, i));
            checkUnderOverFlow(hRsubBC[i]);
        }

        double r[CENTBINS_N], rerr[CENTBINS_N];
        double vobs[CENTBINS_N], vobserr[CENTBINS_N];
        double v[CENTBINS_N], verr[CENTBINS_N];
        for (int i=0; i<CENTBINS_N; i++) {
            cout << "\nCENTRALITY BIN " << i << " : \n";

            double res = GetRes(hRsubAB[i]->GetMean(), hRsubAC[i]->GetMean(), hRsubBC[i]->GetMean());
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
                                  hRsubBC[i]->GetMean(), hRsubBC[i]->GetMeanError());
            vobserr[i] = hVnObs[i]->GetMeanError();
            verr[i] = GetVnError(vobs[i], vobserr[i], r[i], rerr[i]);

            cout << "\n   R" << nCoef << " : " << res << " +- " << rerr[i] << endl;
            cout << "   vobs" << nCoef << " : " << vobs[i] << " +- " << vobserr[i] << endl;
            cout << "   v" << nCoef << " : " << v[i] << " +- " << verr[i] << endl;
        }

        for(int i=0; i<CENTDST_N; i++) {
            if (i==0 || i==1) {
                gR[iFil]->SetPoint(i, double(i+1)*5.0-2.5+iFil, r[i]);
                gR[iFil]->SetPointError(i, 0.0, rerr[i]);

                gVobs[iFil]->SetPoint(i, double(i+1)*5.0-2.5, vobs[i]);
                gVobs[iFil]->SetPointError(i, 0.0, vobserr[i]);

                gV[iFil]->SetPoint(i, double(i+1)*5.0-2.5, v[i]);
                gV[iFil]->SetPointError(i, 0.0, verr[i]);
            } else {
                gR[iFil]->SetPoint(i, double(i)*10.0-5.0+iFil, r[i]);
                gR[iFil]->SetPointError(i, 0.0, rerr[i]);

                gVobs[iFil]->SetPoint(i, double(i)*10.0-5.0, vobs[i]);
                gVobs[iFil]->SetPointError(i, 0.0, vobserr[i]);

                gV[iFil]->SetPoint(i, double(i)*10.0-5.0, v[i]);
                gV[iFil]->SetPointError(i, 0.0, verr[i]);
            }
        }
    }

    //TCanvas *c1 = new TCanvas("c1", "c1", 1200, 500);
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
    //c1->Divide(2, 1);

    //c1->cd(1);
    gR[0]->GetXaxis()->CenterTitle(true);
    gR[0]->GetYaxis()->CenterTitle(true);
    gR[0]->GetYaxis()->SetRangeUser(0.5,1.0);
    gR[0]->SetTitle(Form("Resolution for V0C, granularity on; centrality; R_{%01i}", 2));
    gR[0]->Draw("AP");
    gR[0]->SetMarkerColor(gColor[0]);
    gR[0]->SetMarkerStyle(mMarker);
    gR[0]->SetMarkerSize(mSize);
    gR[0]->Draw("AP");
    c1->Update();

    gR[1]->Draw("P SAME");
    gR[1]->SetMarkerColor(gColor[1]);
    gR[1]->SetMarkerStyle(mMarker);
    gR[1]->SetMarkerSize(mSize);
    gR[1]->Draw("P SAME");
    c1->Update();

    gR[2]->Draw("P SAME");
    gR[2]->SetMarkerColor(gColor[2]);
    gR[2]->SetMarkerStyle(mMarker);
    gR[2]->SetMarkerSize(mSize);
    gR[2]->Draw("P SAME");
    c1->Update();

    TLegend *leg = new TLegend(0.30,0.3,0.55,0.45);
    leg->SetTextSize(0.037);leg->SetBorderSize(0);
    leg->AddEntry(gR[0], "vn scale=1.0", "p");
    leg->AddEntry(gR[1], "vn scale=0.8", "p");
    leg->AddEntry(gR[2], "vn scale=0.65", "p");
    leg->Draw("SAME");

    //c1->cd(2);
    //gV->Draw("AP");
    //gV->GetXaxis()->CenterTitle(true);
    //gV->GetYaxis()->CenterTitle(true);
    //gV->GetYaxis()->SetRangeUser(0.0, 0.25);
    //gV->SetTitle(Form("v_{%01i} for TPC; centrality; v_{%01i}", 2, 2));
    //gV->SetMarkerColor(4);
    //gV->SetMarkerStyle(mrkStyle);
    //gV->SetMarkerSize(mrkSize);
    //gV->Draw("AP");
    //c1->Update();

    //gVobs->Draw("P SAME");
    //gVobs->SetMarkerColor(2);
    //gVobs->GetYaxis()->SetRangeUser(0.0, 0.45);
    //gVobs->SetMarkerStyle(mrkStyle);
    //gVobs->SetMarkerSize(mrkSize);

    //hInput->Draw("HIST SAME");
    //leg->Draw();

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
