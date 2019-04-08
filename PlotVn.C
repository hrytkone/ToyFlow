#include <TMath.h>
#include <TFile.h>

const int nsets = 1;

TString fileNames[nsets] = {"toyFlowGraphs.root"};

TGraphErrors *gVnEP[nsets];
TGraphErrors *gVnEPCorrected[nsets];
TGraphErrors *gVnEPtrue[nsets];
TH1D *hInputFlow[nsets];

TGraphErrors *gR[nsets];
TGraphErrors *gRtrue[nsets];
TGraphErrors *gRcorrected[nsets];

TGraphErrors *gPtBin[nsets];

int mMarker[nsets] = {24};

double mSize = 0.8;

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

        gR[i] = (TGraphErrors*) fIn[i]->Get("gR");
        gRtrue[i] = (TGraphErrors*) fIn[i]->Get("gRtrue");
        gRcorrected[i] = (TGraphErrors*) fIn[i]->Get("gRcorrected");

        gPtBin[i] = (TGraphErrors*) fIn[i]->Get("gPtBin");

        gVnEP[i]->SetMarkerStyle(mMarker[i]);
        gVnEP[i]->SetMarkerColor(i+1);
        gVnEP[i]->SetMarkerSize(mSize);

        gVnEPCorrected[i]->SetMarkerStyle(mMarker[i]);
        gVnEPCorrected[i]->SetMarkerColor(i+2);
        gVnEPCorrected[i]->SetMarkerSize(mSize);

        gVnEPtrue[i]->SetMarkerStyle(mMarker[i]);
        gVnEPtrue[i]->SetMarkerColor(i+3);
        gVnEPtrue[i]->SetMarkerSize(mSize);

        gR[i]->SetMarkerStyle(mMarker[i]+1);
        gR[i]->SetMarkerColor(i+1);
        gR[i]->SetMarkerSize(mSize);

        gRtrue[i]->SetMarkerStyle(mMarker[i]+1);
        gRtrue[i]->SetMarkerColor(i+2);
        gRtrue[i]->SetMarkerSize(mSize);

        gRcorrected[i]->SetMarkerStyle(mMarker[i]+1);
        gRcorrected[i]->SetMarkerColor(i+3);
        gRcorrected[i]->SetMarkerSize(mSize);

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
            gVnEP[i]->Draw("AP");
            leg1->AddEntry(gVnEP[i], "EP", "p");
        }

        gVnEPCorrected[i]->Draw("SAME P");
        leg1->AddEntry(gVnEPCorrected[i], "Corrected EP", "p");

        gVnEPtrue[i]->Draw("SAME P");
        leg1->AddEntry(gVnEPtrue[i], "EP true RP", "p");

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
            gR[i]->Draw("AP");
            leg2->AddEntry(gR[i], "R - sub event method", "p");
        }

        gRtrue[i]->Draw("SAME P");
        leg2->AddEntry(gRtrue[i], "R - true", "p");

        gRcorrected[i]->Draw("SAME P");
        leg2->AddEntry(gRcorrected[i], "R - corrected", "p");
    }

    leg2->Draw("SAME");

    TCanvas *c3 = new TCanvas("c3", "v2(pT)");
    c3->cd();

    TF1 *fVnDist = new TF1("fVnDist", VnDist, 0.0, 6.0, 3);
    fVnDist->SetParameters(2.0, 1.0, 0.15);

    for (i=0; i<nsets; i++) {
        if (i==0) {
            gPtBin[i]->Draw("AP");
        }

        fVnDist->Draw("SAME");
    }

}

double VnDist(double *x, double *p) {
    double pt = x[0];
    double alpha = p[0];
    double beta = p[1];
    double vnMax = p[2];

    double C = vnMax/(TMath::Power(alpha/beta, alpha)*TMath::Exp(-alpha));

    return C*TMath::Power(pt, alpha)*TMath::Exp(-beta*pt);
}
