#include <TMath.h>
#include <TFile.h>

#include "src/JConst.h"

void ErrorExit(TString error, int errnum=1 ){cout<<"ERROR: "<<error<<endl;gSystem->Exit(errnum);}

TString fileName = "toyFlowCentralityGraphs.root";

TGraphErrors *gR[nCoef];
TGraphErrors *gRtrue[nCoef];
TGraphErrors *gVn[nCoef];
TGraphErrors *gVnTrue[nCoef];

TGraphErrors *gRnonuni[nCoef];
TGraphErrors *gRtrueNonuni[nCoef];
TGraphErrors *gVnNonuni[nCoef];
TGraphErrors *gVnTrueNonuni[nCoef];

int mMarker = 24;
double mSize = 1.0;

TFile *fIn;

void PlotCentralityData() {

    fIn = TFile::Open(fileName, "read");

    int i, j;
    for (i=0; i<nCoef; i++) {
        if(fIn==0) ErrorExit(Form("Cannot open file: %s",fileName.Data()));

        gR[i] = (TGraphErrors*) fIn->Get(Form("gRH%02i", i+1));
        gRtrue[i] = (TGraphErrors*) fIn->Get(Form("gRtrueH%02i", i+1));

        gVn[i] = (TGraphErrors*) fIn->Get(Form("gVnH%02i", i+1));
        gVnTrue[i] = (TGraphErrors*) fIn->Get(Form("gVnTrueH%02i", i+1));

        gRnonuni[i] = (TGraphErrors*) fIn->Get(Form("gRnonuniH%02i", i+1));
        gRtrueNonuni[i] = (TGraphErrors*) fIn->Get(Form("gRtrueNonuniH%02i", i+1));

        gVnNonuni[i] = (TGraphErrors*) fIn->Get(Form("gVnNonuniH%02i", i+1));
        gVnTrueNonuni[i] = (TGraphErrors*) fIn->Get(Form("gVnTrueNonuniH%02i", i+1));

        gR[i]->SetMarkerStyle(mMarker);
        gR[i]->SetMarkerColor(1);
        gR[i]->SetMarkerSize(mSize);

        gRtrue[i]->SetMarkerStyle(mMarker+1);
        gRtrue[i]->SetMarkerColor(1);
        gRtrue[i]->SetMarkerSize(mSize);

        gVn[i]->SetMarkerStyle(mMarker);
        gVn[i]->SetMarkerColor(1);
        gVn[i]->SetMarkerSize(mSize);

        gVnTrue[i]->SetMarkerStyle(mMarker+1);
        gVnTrue[i]->SetMarkerColor(1);
        gVnTrue[i]->SetMarkerSize(mSize);

        gRnonuni[i]->SetMarkerStyle(mMarker);
        gRnonuni[i]->SetMarkerColor(2);
        gRnonuni[i]->SetMarkerSize(mSize);

        gRtrueNonuni[i]->SetMarkerStyle(mMarker+1);
        gRtrueNonuni[i]->SetMarkerColor(2);
        gRtrueNonuni[i]->SetMarkerSize(mSize);

        gVnNonuni[i]->SetMarkerStyle(mMarker);
        gVnNonuni[i]->SetMarkerColor(2);
        gVnNonuni[i]->SetMarkerSize(mSize);

        gVnTrueNonuni[i]->SetMarkerStyle(mMarker+1);
        gVnTrueNonuni[i]->SetMarkerColor(2);
        gVnTrueNonuni[i]->SetMarkerSize(mSize);

    }

    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->Divide(2,2);

    for (i=1; i<nCoef; i++) {
        c1->cd(i);
        gVn[i]->Draw("AP");
        gVn[i]->GetYaxis()->SetRangeUser(0.0,0.18);
        gVn[i]->GetXaxis()->CenterTitle(true);
        gVn[i]->GetYaxis()->CenterTitle(true);
        gVn[i]->SetTitle(Form("n=%01i; centrality; v_{n}", i+1));
        gVn[i]->Draw("AP");
        gVnTrue[i]->Draw("SAME P");
        gVnNonuni[i]->Draw("SAME P");
        gVnTrueNonuni[i]->Draw("SAME P");
        c1->Update();
    }

    TLegend *leg1 = new TLegend(0.50,0.6,0.85,0.85,"","brNDC");
    leg1->SetTextSize(0.037);leg1->SetBorderSize(0);

    leg1->AddEntry(gVnTrue[0], "v_{n}, true RP, uniform #phi", "p");
    leg1->AddEntry(gVn[0], "v_{n}, trad. EP, uniform #phi", "p");
    leg1->AddEntry(gVnTrueNonuni[0], "v_{n}, true RP, nonuniform #phi", "p");
    leg1->AddEntry(gVnNonuni[0], "v_{n}, trad. EP, nonuniform #phi", "p");
    leg1->Draw("SAME");
    c1->Draw();

    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->Divide(2,2);

    for (i=1; i<nCoef; i++) {
        c2->cd(i);
        gR[i]->Draw("AP");
        gR[i]->GetYaxis()->SetRangeUser(0.0,1.1);
        gR[i]->GetXaxis()->CenterTitle(true);
        gR[i]->GetYaxis()->CenterTitle(true);
        gR[i]->SetTitle(Form("n=%01i; centrality; R_{n}", i+1));
        gR[i]->Draw("AP");
        gRtrue[i]->Draw("SAME P");
        gRnonuni[i]->Draw("SAME P");
        gRtrueNonuni[i]->Draw("SAME P");
        c2->Update();
    }

    TLegend *leg2 = new TLegend(0.50,0.6,0.85,0.85,"","brNDC");
    leg2->SetTextSize(0.037);leg2->SetBorderSize(0);

    leg2->AddEntry(gRtrue[0], "R true, uniform #phi", "p");
    leg2->AddEntry(gR[0], "R sub event method, uniform #phi", "p");
    leg2->AddEntry(gRtrueNonuni[0], "R true, nonuniform #phi", "p");
    leg2->AddEntry(gRnonuni[0], "R sub event method, nonuniform #phi", "p");
    leg2->Draw("SAME");
    c2->Draw();

}

double VnDist(double *x, double *p) {
    double pt = x[0];
    double alpha = p[0];
    double beta = p[1];
    double vnMax = p[2];
    double C = vnMax/(TMath::Power(alpha/beta, alpha)*TMath::Exp(-alpha));
    return C*TMath::Power(pt, alpha)*TMath::Exp(-beta*pt);
}
