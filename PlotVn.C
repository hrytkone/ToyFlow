#include <TMath.h>
#include <TFile.h>

const int nsets = 1;
const int nCoef = 5;

TString fileNames[nsets] = {"toyFlowGraphs.root"};

TGraphErrors *gVnEP[nsets];
TGraphErrors *gVnEPCorrected[nsets];
TGraphErrors *gVnEPtrue[nsets];
TGraphErrors *gVnEPCorrectedTrue[nsets];
TH1D *hInputFlow[nsets];

TGraphErrors *gR[nsets];
TGraphErrors *gRtrue[nsets];
TGraphErrors *gRcorrected[nsets];
TGraphErrors *gRtrueCorrected[nsets];

TGraphErrors *gVnAlt[nsets];
TGraphErrors *gVnAltCorrected[nsets];

TGraphErrors *gPtBin[nsets];

int mMarker[nsets] = {24};

double mSize = 1.0;

TFile *fIn[nsets];

double VnDist(double *x, double *p);

void PlotVn() {

    int i, j;
    for (i=0; i<nsets; i++) {
        fIn[i] = TFile::Open(fileNames[i], "read");
        gVnEP[i] = (TGraphErrors*) fIn[i]->Get("gVnEP");
        hInputFlow[i] = (TH1D*) fIn[i]->Get("hInputFlow");
        gVnEPtrue[i] = (TGraphErrors*) fIn[i]->Get("gVnEPtrue");
        gVnEPCorrected[i] = (TGraphErrors*) fIn[i]->Get("gVnEPCorrected");
        gVnEPCorrectedTrue[i] = (TGraphErrors*) fIn[i]->Get("gVnEPCorrectedTrue");

        gR[i] = (TGraphErrors*) fIn[i]->Get("gR");
        gRtrue[i] = (TGraphErrors*) fIn[i]->Get("gRtrue");
        gRcorrected[i] = (TGraphErrors*) fIn[i]->Get("gRcorrected");
        gRtrueCorrected[i] = (TGraphErrors*) fIn[i]->Get("gRtrueCorrected");

        gVnAlt[i] = (TGraphErrors*) fIn[i]->Get("gVnAlt");
        gVnAltCorrected[i] = (TGraphErrors*) fIn[i]->Get("gVnAltCorrected");

        gPtBin[i] = (TGraphErrors*) fIn[i]->Get("gPtBin");

        gVnEP[i]->SetMarkerStyle(mMarker[i]);
        gVnEP[i]->SetMarkerColor(i+1);
        gVnEP[i]->SetMarkerSize(mSize);

        gVnEPCorrected[i]->SetMarkerStyle(mMarker[i]);
        gVnEPCorrected[i]->SetMarkerColor(i+2);
        gVnEPCorrected[i]->SetMarkerSize(mSize);

        gVnEPtrue[i]->SetMarkerStyle(mMarker[i]+1);
        gVnEPtrue[i]->SetMarkerColor(i+3);
        gVnEPtrue[i]->SetMarkerSize(mSize);

        gVnEPCorrectedTrue[i]->SetMarkerStyle(mMarker[i]+1);
        gVnEPCorrectedTrue[i]->SetMarkerColor(i+4);
        gVnEPCorrectedTrue[i]->SetMarkerSize(mSize);

        gVnAlt[i]->SetMarkerStyle(mMarker[i]+2);
        gVnAlt[i]->SetMarkerColor(i+6);
        gVnAlt[i]->SetMarkerSize(mSize);

        gVnAltCorrected[i]->SetMarkerStyle(mMarker[i]+2);
        gVnAltCorrected[i]->SetMarkerColor(i+7);
        gVnAltCorrected[i]->SetMarkerSize(mSize);

        gR[i]->SetMarkerStyle(mMarker[i]+1);
        gR[i]->SetMarkerColor(i+1);
        gR[i]->SetMarkerSize(mSize);

        gRcorrected[i]->SetMarkerStyle(mMarker[i]+1);
        gRcorrected[i]->SetMarkerColor(i+2);
        gRcorrected[i]->SetMarkerSize(mSize);

        gRtrue[i]->SetMarkerStyle(mMarker[i]+1);
        gRtrue[i]->SetMarkerColor(i+3);
        gRtrue[i]->SetMarkerSize(mSize);

        gRtrueCorrected[i]->SetMarkerStyle(mMarker[i]+1);
        gRtrueCorrected[i]->SetMarkerColor(i+4);
        gRtrueCorrected[i]->SetMarkerSize(mSize);

        gPtBin[i]->SetMarkerStyle(mMarker[i]);
        gPtBin[i]->SetMarkerColor(i+2);
        gPtBin[i]->SetMarkerSize(mSize);

    }

    TCanvas *c1 = new TCanvas("c1", "vn values with different methods");
    c1->cd();

    TLegend *leg1 = new TLegend(0.65,0.6,0.85,0.75);
    leg1->SetTextSize(0.037);leg1->SetBorderSize(1);

    for (i=0; i<nsets; i++) {
        if (i==0) {
            gVnEPtrue[i]->Draw("AP");
            gVnEPtrue[i]->SetTitle("Vn values with different methods; Vn; n");
            leg1->AddEntry(gVnEPtrue[i], "EP true RP", "p");
        }

        gVnEPCorrectedTrue[i]->Draw("SAME P");
        leg1->AddEntry(gVnEPCorrectedTrue[i], "EP true RP, corrected", "p");

        gVnEP[i]->Draw("SAME P");
        leg1->AddEntry(gVnEP[i], "EP", "p");

        gVnEPCorrected[i]->Draw("SAME P");
        leg1->AddEntry(gVnEPCorrected[i], "EP, corrected", "p");

        gVnAlt[i]->Draw("SAME P");
        leg1->AddEntry(gVnAlt[i], "Alternative EP", "p");

        gVnAltCorrected[i]->Draw("SAME P");
        leg1->AddEntry(gVnAltCorrected[i], "Alternative EP, corrected", "p");

        hInputFlow[i]->Draw("SAME HIST");
        leg1->AddEntry(hInputFlow[i], "Input", "p");
    }

    leg1->Draw("SAME");

    TCanvas *c2 = new TCanvas("c2", "R values with different methods");
    c2->cd();

    TLegend *leg2 = new TLegend(0.65,0.6,0.85,0.75);
    leg2->SetTextSize(0.037);leg2->SetBorderSize(1);

    for (i=0; i<nsets; i++) {
        if (i==0) {
            gRtrue[i]->Draw("AP");
            gRtrue[i]->SetTitle("Resolution parameter; R; n");
            leg2->AddEntry(gRtrue[i], "R true", "p");
        }

        gRtrueCorrected[i]->Draw("SAME P");
        leg2->AddEntry(gRtrueCorrected[i], "R true, corrected", "p");

        gR[i]->Draw("SAME P");
        leg2->AddEntry(gR[i], "R sub event method", "p");

        gRcorrected[i]->Draw("SAME P");
        leg2->AddEntry(gRcorrected[i], "R sub event method, corrected", "p");
    }

    leg2->Draw("SAME");

    TCanvas *c3 = new TCanvas("c3", "v2(pT)");
    c3->cd();

    TF1 *fVnDist = new TF1("fVnDist", VnDist, 0.0, 10.0, 3);
    fVnDist->SetParameters(2.0, 1.0, 0.15);

    for (i=0; i<nsets; i++) {
        if (i==0) {
            gPtBin[i]->Draw("AP");
            gPtBin[i]->SetTitle("v2 as function of pT; v2; pT");
        }

        fVnDist->Draw("SAME");
    }

    TH1D *hPhi = (TH1D*)fIn[0]->Get("hPhi");
    TH1D *hPhiNonuni = (TH1D*)fIn[0]->Get("hPhiNonuni");

    TCanvas *c4 = new TCanvas("c4", "c4");
    c4->Divide(2,1);

    c4->cd(1);
    hPhi->Draw("HIST");

    c4->cd(2);
    hPhiNonuni->Draw("HIST");
}

double VnDist(double *x, double *p) {
    double pt = x[0];
    double alpha = p[0];
    double beta = p[1];
    double vnMax = p[2];
    double C = vnMax/(TMath::Power(alpha/beta, alpha)*TMath::Exp(-alpha));
    return C*TMath::Power(pt, alpha)*TMath::Exp(-beta*pt);
}
