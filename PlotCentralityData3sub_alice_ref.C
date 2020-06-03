#define CENTDST_N 7

#include "src/JConst.h"

void checkUnderOverFlow( TH1 *h );
double GetRes(double rab, double rac, double rbc);
double GetResError(double rab, double rabErr,
                   double rac, double racErr,
                   double rbc, double rbcErr);
double GetVnError(double vobs, double vobsErr, double res, double resErr);

const int nFiles = 4;
const int nFilesRead = 2;
const int nRef = 0;
int gColor[nFiles] = {1,2,4,8};
double scale[nFiles] = {1.0, 0.8, 0.65, 1.0};

int mMarker = 20;
double mSize = 1.0;

TGraphErrors *gR2sub[nFiles];
TGraphErrors *gRtrue2sub[nFiles];
TGraphErrors *gV2sub[nFiles];
TGraphErrors *gVtrue2sub[nFiles];

TGraphErrors *gR[nFiles];
TGraphErrors *gVobs[nFiles];
TGraphErrors *gV[nFiles];

TFile *fIn[nFiles];
TFile *fIn2sub[nFiles];

void PlotCentralityData3sub_alice_ref(TString sFile1 = "testgran0", TString sFile2 = "testgran1", int n = 2, int iVnDet = 6, double multiScale=1.0, double vnScale=1.0) {
    //const double multiScale = 0.8; //0.3, 0.5, 0.8, 1.0
    //const double vnScale = 1.0; //0.65, 0.80, 1.00

    // Näistä luetaan 3-subevent metodiin jutut
    TString fileName[nFilesRead] = {
        //"toyFlow_eta-gap-full-eff-cut_PtDep0_Gran0_Scale0.8/toyFlow_eta-gap-full-eff-cut_PtDep0_Gran0_Scale0.8.root",
        //"toyFlow_eta-gap-full-eff-cut_PtDep0_Gran1_Scale0.8/toyFlow_eta-gap-full-eff-cut_PtDep0_Gran1_Scale0.8.root"
        Form("output/%s/%s.root",sFile1.Data(),sFile1.Data()),
        Form("output/%s/%s.root",sFile2.Data(),sFile2.Data())
    };

    // Nämä luotu MakeCentralityGraphs.C makrolla ja niistä otetaan 2-subevent metodiin asiat
    TString fileName2sub[nFilesRead] = {
        //"toyFlow_eta-gap-full-eff-cut_PtDep0_Gran0_Scale0.8/graphs.root",
        //"toyFlow_eta-gap-full-eff-cut_PtDep0_Gran1_Scale0.8/graphs.root"
        Form("output/%s/%s-plot-output-alice_comp.root",sFile1.Data(),sFile1.Data()),
        Form("output/%s/%s-plot-output-alice_comp.root",sFile2.Data(),sFile2.Data())
    };


    for(int iFil=0; iFil<nFilesRead; iFil++) {
        fIn[iFil] = TFile::Open(fileName[iFil], "read");
        fIn2sub[iFil] = TFile::Open(fileName2sub[iFil], "read");
    }

    Float_t binLower[CENTDST_N+1] = {0, 5, 10, 20, 30, 40, 50, 60};

    double resMeas[CENTDST_N+1] = {0.4862278, 0.6576308, 0.7550672, 0.7798169,
                                   0.7466488, 0.6680212, 0.5434424, 0.383576};

    TH1D *hResMeas = new TH1D("hResMeas", "hResMeas", CENTDST_N, binLower);
    hResMeas->SetLineColor(gColor[0]);
    for (int i=1; i<=CENTDST_N; i++) {
        hResMeas->SetBinContent(i, resMeas[i-1]);
    }

    double centvn[nCoef][CENTDST_N] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0277401, 0.04488324, 0.06521883, 0.08433443, 0.09597485, 0.10087206, 0.09925828},
        //{0.02446622, 0.03929768, 0.05729937, 0.07382241, 0.08401309, 0.08807433, 0.0}, //TPC
        //{0.01974496, 0.03365317, 0.04941064, 0.06331079, 0.07128595, 0.0731899, 0.0}, //V0C
        {0.02039728, 0.02369955, 0.02670301, 0.02950095, 0.03118808, 0.03120636, 0.02918556},
        {0.01013229, 0.01171893, 0.0131265, 0.01479335, 0.0159713, 0.01644628, 0.01535014},
        {0.00415816, 0.00467961, 0.00528238, 0.006501, 0.0068885, 0.00690379, 0.00575251}
    };

    TH1D *hInput[nFilesRead];
    for (int j=0; j<nFilesRead; j++) {
        hInput[j] = new TH1D(Form("hInput%d", j), Form("hInput%d", j), CENTDST_N, binLower);
        hInput[j]->SetLineColor(gColor[0]);
        for (int i=1; i<=CENTDST_N; i++) {
            hInput[j]->SetBinContent(i, vnScale * centvn[n-1][i-1]);
        }
    }

    //const double tpcvn[CENTDST_N] = {0.02446622, 0.03929768, 0.05729937, 0.07382241, 0.08401309, 0.08807433, 0.0};
    //hInput[3] = new TH1D(Form("hInput%d", 3), Form("hInput%d", 3), CENTDST_N, binLower);
    //hInput[3]->SetLineColor(gColor[3]);
    //for (int i=1; i<=CENTDST_N; i++) {
    //    hInput[3]->SetBinContent(i, tpcvn[i-1]);
    //}

    TH1D *hVnObs[CENTBINS_N];
    TH1D *hRsubAB[CENTBINS_N];
    TH1D *hRsubAC[CENTBINS_N];
    TH1D *hRsubBC[CENTBINS_N];

    for(int iFil=0; iFil<nFilesRead; iFil++) {

        // 2SUB METHOD
        gV2sub[iFil] = (TGraphErrors*) fIn2sub[iFil]->Get(Form("gVnH%02i", n));
        gVtrue2sub[iFil] = (TGraphErrors*) fIn2sub[iFil]->Get(Form("gVnTrueH%02i", n));
        gR2sub[iFil] = (TGraphErrors*) fIn2sub[iFil]->Get(Form("gRH%02i", n));
        gRtrue2sub[iFil] = (TGraphErrors*) fIn2sub[iFil]->Get(Form("gRtrueH%02i", n));

        // 3SUB METHOD
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
            if (i==0 || i==1) {
                gR[iFil]->SetPoint(i, double(i+1)*5.0-2.5, r[i]);
                gR[iFil]->SetPointError(i, 0.0, rerr[i]);

                gVobs[iFil]->SetPoint(i, double(i+1)*5.0-2.5, vobs[i]);
                gVobs[iFil]->SetPointError(i, 0.0, vobserr[i]);

                gV[iFil]->SetPoint(i, double(i+1)*5.0-2.5, v[i]);
                gV[iFil]->SetPointError(i, 0.0, verr[i]);
            } else {
                gR[iFil]->SetPoint(i, double(i)*10.0-5.0, r[i]);
                gR[iFil]->SetPointError(i, 0.0, rerr[i]);

                gVobs[iFil]->SetPoint(i, double(i)*10.0-5.0, vobs[i]);
                gVobs[iFil]->SetPointError(i, 0.0, vobserr[i]);

                gV[iFil]->SetPoint(i, double(i)*10.0-5.0, v[i]);
                gV[iFil]->SetPointError(i, 0.0, verr[i]);
            }
        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);

    gRtrue2sub[0]->SetMarkerColor(gColor[0]);
    gRtrue2sub[0]->SetMarkerStyle(mMarker);
    gRtrue2sub[0]->SetMarkerSize(mSize);

    gRtrue2sub[1]->SetMarkerColor(gColor[0]);
    gRtrue2sub[1]->SetMarkerStyle(mMarker+1);
    gRtrue2sub[1]->SetMarkerSize(mSize);

    gR[0]->GetXaxis()->CenterTitle(true);
    gR[0]->GetYaxis()->CenterTitle(true);
    gR[0]->GetYaxis()->SetRangeUser(0.4,1.0);
    gR[0]->GetXaxis()->SetRangeUser(0.0,50.0);
    gR[0]->SetTitle(Form("Resolution of V0C; centrality; R_{%01i}", 2));
    gR[0]->Draw("AP");
    gR[0]->SetMarkerColor(gColor[1]);
    gR[0]->SetFillColor(gColor[1]);
    gR[0]->SetMarkerStyle(mMarker);
    gR[0]->SetMarkerSize(mSize);
    gR[0]->Draw("P SAME");
    c1->Update();

    gR[1]->Draw("P SAME");
    gR[1]->SetMarkerColor(gColor[2]);
    gR[1]->SetFillColor(gColor[2]);
    gR[1]->SetMarkerStyle(mMarker);
    gR[1]->SetMarkerSize(mSize);
    gR[1]->Draw("P SAME");
    c1->Update();

    gR2sub[0]->Draw("P SAME");
    gR2sub[0]->SetMarkerColor(gColor[1]);
    gR2sub[0]->SetMarkerStyle(mMarker+1);
    gR2sub[0]->SetMarkerSize(mSize);
    gR2sub[0]->Draw("P SAME");
    c1->Update();

    gR2sub[1]->Draw("P SAME");
    gR2sub[1]->SetMarkerColor(gColor[2]);
    gR2sub[1]->SetMarkerStyle(mMarker+1);
    gR2sub[1]->SetMarkerSize(mSize);
    gR2sub[1]->Draw("P SAME");
    c1->Update();

    hResMeas->Draw("HIST SAME");
    c1->Update();

    TLegend *leg = new TLegend(0.3,0.25,0.55,0.45,"","brNDC");
    leg->SetTextSize(0.037);leg->SetBorderSize(0);
    leg->SetHeader("#sqrt{s_{NN}} = 5.02 TeV");
    leg->AddEntry(gRtrue2sub[0], "3-sub", "p");
    leg->AddEntry(gRtrue2sub[1], "2-sub", "p");
    leg->AddEntry(gR[0], "Granularity off", "f");
    leg->AddEntry(gR[1], "Granularity on", "f");
    leg->AddEntry(hResMeas, "Measured resolution", "l");
    leg->AddEntry((TH1D*)0x0, Form("Multiplicity scaled by %.2f",multiScale), "l");
    leg->Draw("SAME");
    c1->Update();

    TCanvas *c2 = new TCanvas("c2", "c2", 600, 500);

    gV[0]->Draw("AP");
    gV[0]->GetXaxis()->CenterTitle(true);
    gV[0]->GetYaxis()->CenterTitle(true);
    gV[0]->GetYaxis()->SetRangeUser(0.0, 0.12);
    gV[0]->GetXaxis()->SetRangeUser(0.0, 50.0);
    gV[0]->SetTitle(Form("v_{%01i} from V0C; centrality; v_{%01i}", 2, 2));
    gV[0]->SetMarkerColor(gColor[1]);
    gV[0]->SetFillColor(gColor[1]);
    gV[0]->SetMarkerStyle(mMarker);
    gV[0]->SetMarkerSize(mSize);
    gV[0]->Draw("AP");
    c2->Update();

    gV[1]->Draw("P SAME");
    gV[1]->SetMarkerColor(gColor[2]);
    gV[1]->SetFillColor(gColor[2]);
    gV[1]->SetMarkerStyle(mMarker);
    gV[1]->SetMarkerSize(mSize);
    gV[1]->Draw("P SAME");
    c2->Update();

    gV2sub[0]->Draw("P SAME");
    gV2sub[0]->SetMarkerColor(gColor[1]);
    gV2sub[0]->SetMarkerStyle(mMarker+1);
    gV2sub[0]->SetMarkerSize(mSize);
    gV2sub[0]->Draw("P SAME");
    c2->Update();

    gV2sub[1]->Draw("P SAME");
    gV2sub[1]->SetMarkerColor(gColor[2]);
    gV2sub[1]->SetMarkerStyle(mMarker+1);
    gV2sub[1]->SetMarkerSize(mSize);
    gV2sub[1]->Draw("P SAME");
    c2->Update();

    //hInput[0]->Draw("HIST SAME");
    hInput[1]->Draw("HIST SAME");
    //hInput[2]->Draw("HIST SAME");
    //hInput[3]->Draw("HIST SAME");

    TLegend *leg2 = new TLegend(0.60,0.15,0.85,0.35,"","brNDC");
    leg2->SetTextSize(0.037);leg2->SetBorderSize(0);
    leg2->SetHeader("#sqrt{s_{NN}} = 5.02 TeV");
    leg2->AddEntry(gRtrue2sub[0], "3-sub", "p");
    leg2->AddEntry(gRtrue2sub[1], "2-sub", "p");
    leg2->AddEntry(gV[0], "Granularity off", "p");
    leg2->AddEntry(gV[1], "Granularity on", "p");
    leg2->AddEntry(hInput[0], "Input v_{2}", "l");
    leg2->AddEntry(hInput[1], Form("v_{%d} input, scale=%.2f", 2, vnScale), "l");
    leg2->AddEntry((TH1D*)0x0, Form("Multiplicity scaled by %.2f",multiScale), "l");
    leg2->Draw("SAME");
    c2->Update();

    c1->SaveAs(Form("output/%s/%s-reso.pdf",sFile1.Data(),sFile1.Data()));
    c2->SaveAs(Form("output/%s/%s-v2.pdf",sFile1.Data(),sFile1.Data()));
    c1->SaveAs(Form("output/%s/%s-reso.pdf",sFile2.Data(),sFile2.Data()));
    c2->SaveAs(Form("output/%s/%s-v2.pdf",sFile2.Data(),sFile2.Data()));
}

//______________________________________________________________________________

void checkUnderOverFlow( TH1 *h ) {
        if(h->GetBinContent(0)>0) cout << h->GetName() << " underflow bin not empty: " << h->GetBinContent(0) << endl;
        if(h->GetBinContent(h->GetXaxis()->GetNbins()+1)>0) cout << h->GetName() << " overflow bin not empty: " << h->GetBinContent(h->GetXaxis()->GetNbins()+1) << endl;
}

double GetRes(double rab, double rac, double rbc) {
    return TMath::Sqrt((rab * rac)/rbc);
}

double GetResError(double rab, double rabErr,
                   double rac, double racErr,
                   double rbc, double rbcErr) {
    return 0.5*GetRes(rab, rac, rbc)*TMath::Sqrt(TMath::Power(rabErr/rab, 2.0) +
    TMath::Power(racErr/rac, 2.0) + TMath::Power(rbcErr/rbc, 2.0));
}

double GetVnError(double vobs, double vobsErr, double res, double resErr) {
    return (vobs/res)*TMath::Sqrt(TMath::Power(vobsErr/vobs, 2.0) +
    TMath::Power(resErr/res, 2.0));
}
