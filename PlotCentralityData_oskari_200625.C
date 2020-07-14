#define CENTDST_N 7

#include "src/JConst.h"

void checkUnderOverFlow( TH1 *h );
double GetRes(double rab, double rac, double rbc);
double GetResError(double rab, double rabErr,
                   double rac, double racErr,
                   double rbc, double rbcErr);
double GetVnError(double vobs, double vobsErr, double res, double resErr);

const int nFiles = 14;
const int nFilesRead = 14;
const int nRef = 0;
int gColor[5] = {1,2,4,8,9};
//double scale[4] = {1.0, 0.8, 0.65, 1.0};

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

void PlotCentralityData_oskari_200625(int n = 2, int iVnDet = 6, double multiScale=1.0, double vnScale=1.0) {
    //const double multiScale = 0.8; //0.3, 0.5, 0.8, 1.0
    //const double vnScale = 1.0; //0.65, 0.80, 1.00

    // Näistä luetaan 3-subevent metodiin jutut
    TString fileName[nFilesRead] = {
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0.root", // Ideal
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran1_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran1_Scale1.00_ptCut0.root", // Granularity on
        "output/toyFlow_20200630_0.8Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_0.8Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0.root", // Multi -20 %
        "output/toyFlow_20200630_0.5Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_0.5Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0.root", // Multi -50 %
        "output/toyFlow_20200630_0.3Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_0.3Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0.root", // Multi -70 %
        "output/toyFlow_20200630_1.0Multi_0.1extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.1extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0.root", // Conversions 10 %
        "output/toyFlow_20200630_1.0Multi_0.2extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.2extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0.root", // Conversions 20 %
        "output/toyFlow_20200630_1.0Multi_0.5extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.5extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0.root", // Conversions 50 %
        "output/toyFlow_20200630_1.0Multi_0.8extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.8extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0.root", // Conversions 80 %
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.1bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.1bDecay_PtDep0_Gran0_Scale1.00_ptCut0.root", // Decays 10 %
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.2bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.2bDecay_PtDep0_Gran0_Scale1.00_ptCut0.root", // Decays 20 %
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale0.80_ptCut0/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale0.80_ptCut0.root", // v_n -20 %
        "output/toyFlow_20200630_0.8Multi_0.2extraConvPart_0.0bDecay_PtDep0_Gran1_Scale0.80_ptCut0/toyFlow_20200630_0.8Multi_0.2extraConvPart_0.0bDecay_PtDep0_Gran1_Scale0.80_ptCut0.root", // Best guess
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut1/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut1.root"  // Pt cut for TPC
    };

    // Nämä luotu MakeCentralityGraphs.C makrolla ja niistä otetaan 2-subevent metodiin asiat
    TString fileName2sub[nFilesRead] = {
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0-plot-output-alice_comp.root", // Ideal
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran1_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran1_Scale1.00_ptCut0-plot-output-alice_comp.root", // Granularity on
        "output/toyFlow_20200630_0.8Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_0.8Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0-plot-output-alice_comp.root", // Multi -20 %
        "output/toyFlow_20200630_0.5Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_0.5Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0-plot-output-alice_comp.root", // Multi -50 %
        "output/toyFlow_20200630_0.3Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_0.3Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0-plot-output-alice_comp.root", // Multi -70 %
        "output/toyFlow_20200630_1.0Multi_0.1extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.1extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0-plot-output-alice_comp.root", // Conversions 10 %
        "output/toyFlow_20200630_1.0Multi_0.2extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.2extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0-plot-output-alice_comp.root", // Conversions 20 %
        "output/toyFlow_20200630_1.0Multi_0.5extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.5extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0-plot-output-alice_comp.root", // Conversions 50 %
        "output/toyFlow_20200630_1.0Multi_0.8extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.8extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut0-plot-output-alice_comp.root", // Conversions 80 %
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.1bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.1bDecay_PtDep0_Gran0_Scale1.00_ptCut0-plot-output-alice_comp.root", // Decays 10 %
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.2bDecay_PtDep0_Gran0_Scale1.00_ptCut0/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.2bDecay_PtDep0_Gran0_Scale1.00_ptCut0-plot-output-alice_comp.root", // Decays 20 %
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale0.80_ptCut0/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale0.80_ptCut0-plot-output-alice_comp.root", // v_n -20 %
        "output/toyFlow_20200630_0.8Multi_0.2extraConvPart_0.0bDecay_PtDep0_Gran1_Scale0.80_ptCut0/toyFlow_20200630_0.8Multi_0.2extraConvPart_0.0bDecay_PtDep0_Gran1_Scale0.80_ptCut0-plot-output-alice_comp.root", // Best guess
        "output/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut1/toyFlow_20200630_1.0Multi_0.0extraConvPart_0.0bDecay_PtDep0_Gran0_Scale1.00_ptCut1-plot-output-alice_comp.root"  // Pt cut for TPC
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

    gR2sub[0]->Draw("P SAME");
    gR2sub[0]->SetMarkerColor(gColor[1]);
    gR2sub[0]->SetMarkerStyle(mMarker+1);
    gR2sub[0]->SetMarkerSize(mSize);
    gR2sub[0]->Draw("P SAME");
    c1->Update();

    hResMeas->Draw("HIST SAME");
    c1->Update();

    TLegend *leg = new TLegend(0.3,0.15,0.55,0.40,"Ideal","brNDC");
    leg->SetTextSize(0.037);leg->SetBorderSize(0);
    leg->SetHeader("Ideal, #sqrt{s_{NN}} = 5.02 TeV");
    leg->AddEntry(gR[0], "3-sub", "p");
    leg->AddEntry(gR2sub[0], "2-sub", "p");
    leg->AddEntry(hResMeas, "Measured resolution", "l");
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

    gV2sub[0]->Draw("P SAME");
    gV2sub[0]->SetMarkerColor(gColor[1]);
    gV2sub[0]->SetMarkerStyle(mMarker+1);
    gV2sub[0]->SetMarkerSize(mSize);
    gV2sub[0]->Draw("P SAME");
    c2->Update();

    //hInput[0]->Draw("HIST SAME");
    hInput[1]->Draw("HIST SAME");
    //hInput[2]->Draw("HIST SAME");
    //hInput[3]->Draw("HIST SAME");

    TLegend *leg2 = new TLegend(0.40,0.15,0.65,0.35,"","brNDC");
    leg2->SetTextSize(0.037);leg2->SetBorderSize(0);
    leg2->SetHeader("Ideal, #sqrt{s_{NN}} = 5.02 TeV");
    leg2->AddEntry(gV[0], "3-sub", "p");
    leg2->AddEntry(gV2sub[0], "2-sub", "p");
    leg2->AddEntry(hInput[0], "Input v_{2}", "l");
    leg2->Draw("SAME");
    c2->Update();



    TCanvas *cgran = new TCanvas("cgran", "cgran", 600, 500);

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
    cgran->Update();

    gR[1]->Draw("P SAME");
    gR[1]->SetMarkerColor(gColor[2]);
    gR[1]->SetFillColor(gColor[2]);
    gR[1]->SetMarkerStyle(mMarker);
    gR[1]->SetMarkerSize(mSize);
    gR[1]->Draw("P SAME");
    cgran->Update();

    gR2sub[0]->Draw("P SAME");
    gR2sub[0]->SetMarkerColor(gColor[1]);
    gR2sub[0]->SetMarkerStyle(mMarker+1);
    gR2sub[0]->SetMarkerSize(mSize);
    gR2sub[0]->Draw("P SAME");
    cgran->Update();

    gR2sub[1]->Draw("P SAME");
    gR2sub[1]->SetMarkerColor(gColor[2]);
    gR2sub[1]->SetMarkerStyle(mMarker+1);
    gR2sub[1]->SetMarkerSize(mSize);
    gR2sub[1]->Draw("P SAME");
    cgran->Update();

    hResMeas->Draw("HIST SAME");
    cgran->Update();

    TLegend *leg3 = new TLegend(0.3,0.15,0.55,0.40,"","brNDC");
    leg3->SetTextSize(0.037);leg3->SetBorderSize(0);
    leg3->SetHeader("Granularity, #sqrt{s_{NN}} = 5.02 TeV");
    leg3->AddEntry(gRtrue2sub[0], "3-sub", "p");
    leg3->AddEntry(gRtrue2sub[1], "2-sub", "p");
    leg3->AddEntry(gR[0], "Granularity off", "f");
    leg3->AddEntry(gR[1], "Granularity on", "f");
    leg3->AddEntry(hResMeas, "Measured resolution", "l");
    leg3->Draw("SAME");
    cgran->Update();


    TCanvas *cgran2 = new TCanvas("cgran2", "cgran2", 600, 500);

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
    cgran2->Update();

    gV[1]->Draw("P SAME");
    gV[1]->SetMarkerColor(gColor[2]);
    gV[1]->SetFillColor(gColor[2]);
    gV[1]->SetMarkerStyle(mMarker);
    gV[1]->SetMarkerSize(mSize);
    gV[1]->Draw("P SAME");
    cgran2->Update();

    gV2sub[0]->Draw("P SAME");
    gV2sub[0]->SetMarkerColor(gColor[1]);
    gV2sub[0]->SetMarkerStyle(mMarker+1);
    gV2sub[0]->SetMarkerSize(mSize);
    gV2sub[0]->Draw("P SAME");
    cgran2->Update();

    gV2sub[1]->Draw("P SAME");
    gV2sub[1]->SetMarkerColor(gColor[2]);
    gV2sub[1]->SetMarkerStyle(mMarker+1);
    gV2sub[1]->SetMarkerSize(mSize);
    gV2sub[1]->Draw("P SAME");
    cgran2->Update();

    //hInput[0]->Draw("HIST SAME");
    hInput[1]->Draw("HIST SAME");
    //hInput[2]->Draw("HIST SAME");
    //hInput[3]->Draw("HIST SAME");

    TLegend *leg4 = new TLegend(0.40,0.15,0.65,0.35,"","brNDC");
    leg4->SetTextSize(0.037);leg4->SetBorderSize(0);
    leg4->SetHeader("Granularity, #sqrt{s_{NN}} = 5.02 TeV");
    leg4->AddEntry(gRtrue2sub[0], "3-sub", "p");
    leg4->AddEntry(gRtrue2sub[1], "2-sub", "p");
    leg4->AddEntry(gV[0], "Granularity off", "p");
    leg4->AddEntry(gV[1], "Granularity on", "p");
    leg4->AddEntry(hInput[0], "Input v_{2}", "l");
    leg4->Draw("SAME");
    cgran2->Update();



    TCanvas *cMultiRed = new TCanvas("cMultiRed", "cMultiRed", 600, 500);

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
    cMultiRed->Update();

    gR[2]->Draw("P SAME");
    gR[2]->SetMarkerColor(gColor[1]);
    gR[2]->SetMarkerStyle(mMarker+1);
    gR[2]->SetMarkerSize(mSize);
    gR[2]->Draw("P SAME");
    cMultiRed->Update();

    gR[3]->Draw("P SAME");
    gR[3]->SetMarkerColor(gColor[2]);
    gR[3]->SetMarkerStyle(mMarker+2);
    gR[3]->SetMarkerSize(mSize);
    gR[3]->Draw("P SAME");
    cMultiRed->Update();

    gR[4]->Draw("P SAME");
    gR[4]->SetMarkerColor(gColor[3]);
    gR[4]->SetMarkerStyle(mMarker+3);
    gR[4]->SetMarkerSize(mSize);
    gR[4]->Draw("P SAME");
    cMultiRed->Update();

    hResMeas->Draw("HIST SAME");
    cMultiRed->Update();

    TLegend *legMultiRed = new TLegend(0.3,0.15,0.55,0.40,"Ideal","brNDC");
    legMultiRed->SetTextSize(0.037);legMultiRed->SetBorderSize(0);
    legMultiRed->SetHeader("Multiplicity reduction, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legMultiRed->AddEntry(gR[0], "Ideal", "p");
    legMultiRed->AddEntry(gR[2], "Multiplicity reduced 20 %", "p");
    legMultiRed->AddEntry(gR[3], "Multiplicity reduced 50 %", "p");
    legMultiRed->AddEntry(gR[4], "Multiplicity reduced 80 %", "p");
    legMultiRed->AddEntry(hResMeas, "Measured resolution", "l");
    legMultiRed->Draw("SAME");
    cMultiRed->Update();



    TCanvas *cMultiRed2 = new TCanvas("cMultiRed2", "cMultiRed2", 600, 500);

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
    cMultiRed2->Update();

    gV[2]->Draw("P SAME");
    gV[2]->SetMarkerColor(gColor[1]);
    gV[2]->SetMarkerStyle(mMarker+1);
    gV[2]->SetMarkerSize(mSize);
    gV[2]->Draw("P SAME");
    cMultiRed2->Update();

    gV[3]->Draw("P SAME");
    gV[3]->SetMarkerColor(gColor[2]);
    gV[3]->SetMarkerStyle(mMarker+2);
    gV[3]->SetMarkerSize(mSize);
    gV[3]->Draw("P SAME");
    cMultiRed2->Update();

    gV[4]->Draw("P SAME");
    gV[4]->SetMarkerColor(gColor[3]);
    gV[4]->SetMarkerStyle(mMarker+3);
    gV[4]->SetMarkerSize(mSize);
    gV[4]->Draw("P SAME");
    cMultiRed2->Update();

    //hInput[0]->Draw("HIST SAME");
    hInput[1]->Draw("HIST SAME");
    //hInput[2]->Draw("HIST SAME");
    //hInput[3]->Draw("HIST SAME");

    TLegend *legMultiRed2 = new TLegend(0.25,0.15,0.50,0.35,"","brNDC");
    legMultiRed2->SetTextSize(0.037);legMultiRed2->SetBorderSize(0);
    legMultiRed2->SetHeader("Multiplicity reduction, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legMultiRed2->AddEntry(gV[0], "Ideal", "p");
    legMultiRed2->AddEntry(gV[2], "Multiplicity reduced 20 %", "p");
    legMultiRed2->AddEntry(gV[3], "Multiplicity reduced 50 %", "p");
    legMultiRed2->AddEntry(gV[4], "Multiplicity reduced 80 %", "p");
    legMultiRed2->AddEntry(hInput[0], "Input v_{2}", "l");
    legMultiRed2->Draw("SAME");
    cMultiRed2->Update();



    TCanvas *cSecondaries = new TCanvas("cSecondaries", "cSecondaries", 600, 500);

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
    cSecondaries->Update();

    gR[5]->Draw("P SAME");
    gR[5]->SetMarkerColor(gColor[1]);
    gR[5]->SetMarkerStyle(mMarker+1);
    gR[5]->SetMarkerSize(mSize);
    gR[5]->Draw("P SAME");
    cSecondaries->Update();

    gR[6]->Draw("P SAME");
    gR[6]->SetMarkerColor(gColor[2]);
    gR[6]->SetMarkerStyle(mMarker+2);
    gR[6]->SetMarkerSize(mSize);
    gR[6]->Draw("P SAME");
    cSecondaries->Update();

    gR[7]->Draw("P SAME");
    gR[7]->SetMarkerColor(gColor[3]);
    gR[7]->SetMarkerStyle(mMarker+3);
    gR[7]->SetMarkerSize(mSize);
    gR[7]->Draw("P SAME");
    cSecondaries->Update();

    gR[8]->Draw("P SAME");
    gR[8]->SetMarkerColor(gColor[4]);
    gR[8]->SetMarkerStyle(mMarker+4);
    gR[8]->SetMarkerSize(mSize);
    gR[8]->Draw("P SAME");
    cSecondaries->Update();

    hResMeas->Draw("HIST SAME");
    cSecondaries->Update();

    TLegend *legSecondaries = new TLegend(0.3,0.15,0.45,0.35,"Ideal","brNDC");
    legSecondaries->SetTextSize(0.037);legSecondaries->SetBorderSize(0);
    legSecondaries->SetHeader("Secondaries, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legSecondaries->AddEntry(gR[0], "Ideal", "p");
    legSecondaries->AddEntry(gR[5], "Secondaries added 10 %", "p");
    legSecondaries->AddEntry(gR[6], "Secondaries added 20 %", "p");
    legSecondaries->AddEntry(gR[7], "Secondaries added 50 %", "p");
    legSecondaries->AddEntry(gR[8], "Secondaries added 80 %", "p");
    legSecondaries->AddEntry(hResMeas, "Measured resolution", "l");
    legSecondaries->Draw("SAME");
    cSecondaries->Update();



    TCanvas *cSecondaries2 = new TCanvas("cSecondaries2", "cSecondaries2", 600, 500);

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
    cSecondaries2->Update();

    gV[5]->Draw("P SAME");
    gV[5]->SetMarkerColor(gColor[1]);
    gV[5]->SetMarkerStyle(mMarker+1);
    gV[5]->SetMarkerSize(mSize);
    gV[5]->Draw("P SAME");
    cSecondaries2->Update();

    gV[6]->Draw("P SAME");
    gV[6]->SetMarkerColor(gColor[2]);
    gV[6]->SetMarkerStyle(mMarker+2);
    gV[6]->SetMarkerSize(mSize);
    gV[6]->Draw("P SAME");
    cSecondaries2->Update();

    gV[7]->Draw("P SAME");
    gV[7]->SetMarkerColor(gColor[3]);
    gV[7]->SetMarkerStyle(mMarker+3);
    gV[7]->SetMarkerSize(mSize);
    gV[7]->Draw("P SAME");
    cSecondaries2->Update();

    gV[8]->Draw("P SAME");
    gV[8]->SetMarkerColor(gColor[4]);
    gV[8]->SetMarkerStyle(mMarker+4);
    gV[8]->SetMarkerSize(mSize);
    gV[8]->Draw("P SAME");
    cSecondaries2->Update();

    //hInput[0]->Draw("HIST SAME");
    hInput[1]->Draw("HIST SAME");
    //hInput[2]->Draw("HIST SAME");
    //hInput[3]->Draw("HIST SAME");

    TLegend *legSecondaries2 = new TLegend(0.15,0.65,0.40,0.85,"","brNDC");
    legSecondaries2->SetTextSize(0.037);legSecondaries2->SetBorderSize(0);
    legSecondaries2->SetHeader("Secondaries, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legSecondaries2->AddEntry(gV[0], "Ideal", "p");
    legSecondaries2->AddEntry(gV[5], "Secondaries added 10 %", "p");
    legSecondaries2->AddEntry(gV[6], "Secondaries added 20 %", "p");
    legSecondaries2->AddEntry(gV[7], "Secondaries added 50 %", "p");
    legSecondaries2->AddEntry(gV[8], "Secondaries added 80 %", "p");
    legSecondaries2->AddEntry(hInput[0], "Input v_{2}", "l");
    legSecondaries2->Draw("SAME");
    cSecondaries2->Update();



    TCanvas *cDecays = new TCanvas("cDecays", "cDecays", 600, 500);

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
    cDecays->Update();

    gR[9]->Draw("P SAME");
    gR[9]->SetMarkerColor(gColor[1]);
    gR[9]->SetMarkerStyle(mMarker+1);
    gR[9]->SetMarkerSize(mSize);
    gR[9]->Draw("P SAME");
    cDecays->Update();

    gR[10]->Draw("P SAME");
    gR[10]->SetMarkerColor(gColor[2]);
    gR[10]->SetMarkerStyle(mMarker+2);
    gR[10]->SetMarkerSize(mSize);
    gR[10]->Draw("P SAME");
    cDecays->Update();

    hResMeas->Draw("HIST SAME");
    cDecays->Update();

    TLegend *legDecays = new TLegend(0.3,0.15,0.55,0.40,"Ideal","brNDC");
    legDecays->SetTextSize(0.037);legDecays->SetBorderSize(0);
    legDecays->SetHeader("Decays, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legDecays->AddEntry(gR[0], "Ideal", "p");
    legDecays->AddEntry(gR[9], "10 % of all particles decay", "p");
    legDecays->AddEntry(gR[10],"20 % of all particles decay", "p");
    legDecays->AddEntry(hResMeas, "Measured resolution", "l");
    legDecays->Draw("SAME");
    cDecays->Update();



    TCanvas *cDecays2 = new TCanvas("cDecays2", "cDecays2", 600, 500);

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
    cDecays2->Update();

    gV[9]->Draw("P SAME");
    gV[9]->SetMarkerColor(gColor[1]);
    gV[9]->SetMarkerStyle(mMarker+1);
    gV[9]->SetMarkerSize(mSize);
    gV[9]->Draw("P SAME");
    cDecays2->Update();

    gV[10]->Draw("P SAME");
    gV[10]->SetMarkerColor(gColor[2]);
    gV[10]->SetMarkerStyle(mMarker+2);
    gV[10]->SetMarkerSize(mSize);
    gV[10]->Draw("P SAME");
    cDecays2->Update();

    //hInput[0]->Draw("HIST SAME");
    hInput[1]->Draw("HIST SAME");
    //hInput[2]->Draw("HIST SAME");
    //hInput[3]->Draw("HIST SAME");

    TLegend *legDecays2 = new TLegend(0.40,0.15,0.65,0.35,"","brNDC");
    legDecays2->SetTextSize(0.037);legDecays2->SetBorderSize(0);
    legDecays2->SetHeader("Decays, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legDecays2->AddEntry(gV[0], "Ideal", "p");
    legDecays2->AddEntry(gV[9], "10 % of all particles decay", "p");
    legDecays2->AddEntry(gV[10],"20 % of all particles decay", "p");
    legDecays2->AddEntry(hInput[0], "Input v_{2}", "l");
    legDecays2->Draw("SAME");
    cDecays2->Update();



    TCanvas *cVnScale = new TCanvas("cVnScale", "cVnScale", 600, 500);

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
    cVnScale->Update();

    gR[11]->Draw("P SAME");
    gR[11]->SetMarkerColor(gColor[1]);
    gR[11]->SetMarkerStyle(mMarker+1);
    gR[11]->SetMarkerSize(mSize);
    gR[11]->Draw("P SAME");
    cVnScale->Update();

    hResMeas->Draw("HIST SAME");
    cVnScale->Update();

    TLegend *legVnScale = new TLegend(0.3,0.15,0.55,0.40,"Ideal","brNDC");
    legVnScale->SetTextSize(0.037);legVnScale->SetBorderSize(0);
    legVnScale->SetHeader("Vn scaled, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legVnScale->AddEntry(gR[0], "Ideal", "p");
    legVnScale->AddEntry(gR[11], "Vn scaled down 20 %", "p");
    legVnScale->AddEntry(hResMeas, "Measured resolution", "l");
    legVnScale->Draw("SAME");
    cVnScale->Update();



    TCanvas *cVnScale2 = new TCanvas("cVnScale2", "cVnScale2", 600, 500);

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
    cVnScale2->Update();

    gV[11]->Draw("P SAME");
    gV[11]->SetMarkerColor(gColor[1]);
    gV[11]->SetMarkerStyle(mMarker+1);
    gV[11]->SetMarkerSize(mSize);
    gV[11]->Draw("P SAME");
    cVnScale2->Update();

    //hInput[0]->Draw("HIST SAME");
    hInput[1]->Draw("HIST SAME");
    //hInput[2]->Draw("HIST SAME");
    //hInput[3]->Draw("HIST SAME");

    TLegend *legVnScale2 = new TLegend(0.40,0.15,0.65,0.35,"","brNDC");
    legVnScale2->SetTextSize(0.037);legVnScale2->SetBorderSize(0);
    legVnScale2->SetHeader("Vn scaled, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legVnScale2->AddEntry(gV[0], "Ideal", "p");
    legVnScale2->AddEntry(gV[11], "Vn scaled down 20 %", "p");
    legVnScale2->AddEntry(hInput[0], "Input v_{2}", "l");
    legVnScale2->Draw("SAME");
    cVnScale2->Update();




    TCanvas *cFinal = new TCanvas("cFinal", "cFinal", 600, 500);

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
    cFinal->Update();

    gR[12]->Draw("P SAME");
    gR[12]->SetMarkerColor(gColor[1]);
    gR[12]->SetMarkerStyle(mMarker+1);
    gR[12]->SetMarkerSize(mSize);
    gR[12]->Draw("P SAME");
    cFinal->Update();

    hResMeas->Draw("HIST SAME");
    cFinal->Update();

    TLegend *legFinal = new TLegend(0.3,0.15,0.55,0.40,"Ideal","brNDC");
    legFinal->SetTextSize(0.037);legFinal->SetBorderSize(0);
    legFinal->SetHeader("Best estimate, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legFinal->AddEntry(gR[0], "Ideal", "p");
    legFinal->AddEntry(gR[12], "tracking -20%, secondaries +20%", "p");
    legFinal->AddEntry((TH1D*)0x0, "vn scaled -20%, granularity on", "");
    legFinal->AddEntry(hResMeas, "Measured resolution", "l");
    legFinal->Draw("SAME");
    cFinal->Update();



    TCanvas *cFinal2 = new TCanvas("cFinal2", "cFinal2", 600, 500);

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
    cFinal2->Update();

    gV[12]->Draw("P SAME");
    gV[12]->SetMarkerColor(gColor[1]);
    gV[12]->SetMarkerStyle(mMarker+1);
    gV[12]->SetMarkerSize(mSize);
    gV[12]->Draw("P SAME");
    cFinal2->Update();

    //hInput[0]->Draw("HIST SAME");
    hInput[1]->Draw("HIST SAME");
    //hInput[2]->Draw("HIST SAME");
    //hInput[3]->Draw("HIST SAME");

    TLegend *legFinal2 = new TLegend(0.40,0.15,0.65,0.35,"","brNDC");
    legFinal2->SetTextSize(0.037);legFinal2->SetBorderSize(0);
    legFinal2->SetHeader("Best estimate, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legFinal2->AddEntry(gV[0], "Ideal", "p");
    legFinal2->AddEntry(gV[12], "tracking -20%, secondaries +20%", "p");
    legFinal2->AddEntry((TH1D*)0x0, "vn scaled -20%, granularity on", "");
    legFinal2->AddEntry(hInput[0], "Input v_{2}", "l");
    legFinal2->Draw("SAME");
    cFinal2->Update();




    TCanvas *cPtCut = new TCanvas("cPtCut", "cPtCut", 600, 500);

    gR2sub[0]->GetXaxis()->CenterTitle(true);
    gR2sub[0]->GetYaxis()->CenterTitle(true);
    gR2sub[0]->GetYaxis()->SetRangeUser(0.4,1.0);
    gR2sub[0]->GetXaxis()->SetRangeUser(0.0,50.0);
    gR2sub[0]->SetTitle(Form("Resolution of V0C; centrality; R_{%01i}", 2));
    gR2sub[0]->Draw("AP");
    gR2sub[0]->SetMarkerColor(gColor[1]);
    gR2sub[0]->SetFillColor(gColor[1]);
    gR2sub[0]->SetMarkerStyle(mMarker);
    gR2sub[0]->SetMarkerSize(mSize);
    gR2sub[0]->Draw("P SAME");
    cPtCut->Update();

    gR2sub[13]->Draw("P SAME");
    gR2sub[13]->SetMarkerColor(gColor[1]);
    gR2sub[13]->SetMarkerStyle(mMarker+1);
    gR2sub[13]->SetMarkerSize(mSize);
    gR2sub[13]->Draw("P SAME");
    cPtCut->Update();

    hResMeas->Draw("HIST SAME");
    cPtCut->Update();

    TLegend *legPtCut = new TLegend(0.3,0.15,0.55,0.40,"Ideal","brNDC");
    legPtCut->SetTextSize(0.037);legPtCut->SetBorderSize(0);
    legPtCut->SetHeader("Pt cut, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legPtCut->AddEntry(gR2sub[0], "Ideal", "p");
    legPtCut->AddEntry(gR2sub[13], "Pt cut < 150 MeV in TPC", "p");
    legPtCut->AddEntry(hResMeas, "Measured resolution", "l");
    legPtCut->Draw("SAME");
    cPtCut->Update();



    TCanvas *cPtCut2 = new TCanvas("cPtCut2", "cPtCut2", 600, 500);

    gV2sub[0]->Draw("AP");
    gV2sub[0]->GetXaxis()->CenterTitle(true);
    gV2sub[0]->GetYaxis()->CenterTitle(true);
    gV2sub[0]->GetYaxis()->SetRangeUser(0.0, 0.12);
    gV2sub[0]->GetXaxis()->SetRangeUser(0.0, 50.0);
    gV2sub[0]->SetTitle(Form("v_{%01i} from V0C; centrality; v_{%01i}", 2, 2));
    gV2sub[0]->SetMarkerColor(gColor[1]);
    gV2sub[0]->SetFillColor(gColor[1]);
    gV2sub[0]->SetMarkerStyle(mMarker);
    gV2sub[0]->SetMarkerSize(mSize);
    gV2sub[0]->Draw("AP");
    cPtCut2->Update();

    gV2sub[13]->Draw("P SAME");
    gV2sub[13]->SetMarkerColor(gColor[1]);
    gV2sub[13]->SetMarkerStyle(mMarker+1);
    gV2sub[13]->SetMarkerSize(mSize);
    gV2sub[13]->Draw("P SAME");
    cPtCut2->Update();

    //hInput[0]->Draw("HIST SAME");
    hInput[1]->Draw("HIST SAME");
    //hInput[2]->Draw("HIST SAME");
    //hInput[3]->Draw("HIST SAME");

    TLegend *legPtCut2 = new TLegend(0.40,0.15,0.65,0.35,"","brNDC");
    legPtCut2->SetTextSize(0.037);legPtCut2->SetBorderSize(0);
    legPtCut2->SetHeader("Pt cut, 3-subevt, #sqrt{s_{NN}} = 5.02 TeV");
    legPtCut2->AddEntry(gV2sub[0], "Ideal", "p");
    legPtCut2->AddEntry(gV2sub[13], "Pt cut < 150 MeV in TPC", "p");
    legPtCut2->AddEntry(hInput[0], "Input v_{2}", "l");
    legPtCut2->Draw("SAME");
    cPtCut2->Update();



    c1->SaveAs(Form("figures/oskari_200625/Ideal_R2.pdf"));
    c2->SaveAs(Form("figures/oskari_200625/Ideal_v2.pdf"));
    cgran->SaveAs(Form("figures/oskari_200625/Gran_R2.pdf"));
    cgran2->SaveAs(Form("figures/oskari_200625/Gran_v2.pdf"));
    cMultiRed->SaveAs(Form("figures/oskari_200625/Ineff_R2.pdf"));
    cMultiRed2->SaveAs(Form("figures/oskari_200625/Ineff_v2.pdf"));
    cSecondaries->SaveAs(Form("figures/oskari_200625/Second_R2.pdf"));
    cSecondaries2->SaveAs(Form("figures/oskari_200625/Second_v2.pdf"));
    cDecays->SaveAs(Form("figures/oskari_200625/Decays_R2.pdf"));
    cDecays2->SaveAs(Form("figures/oskari_200625/Decays_v2.pdf"));
    cVnScale->SaveAs(Form("figures/oskari_200625/VnScale_R2.pdf"));
    cVnScale2->SaveAs(Form("figures/oskari_200625/VnScale_v2.pdf"));
    cFinal->SaveAs(Form("figures/oskari_200625/Final-R2.pdf"));
    cFinal2->SaveAs(Form("figures/oskari_200625/Final-v2.pdf"));
    cPtCut->SaveAs(Form("figures/oskari_200625/PtCut-R2.pdf"));
    cPtCut2->SaveAs(Form("figures/oskari_200625/PtCut-v2.pdf"));
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
