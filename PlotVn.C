#include <TMath.h>
#include <TFile.h>

const int nsets = 1;

void ErrorExit(TString error, int errnum=1 ){cout<<"ERROR: "<<error<<endl;gSystem->Exit(errnum);}

TString fileNames[nsets] = {"toyFlowGraphs.root"};

TH1D *hInputFlow[nsets];

TGraphErrors *gR[nsets];
TGraphErrors *gRtrue[nsets];
TGraphErrors *gVn[nsets];
TGraphErrors *gVnTrue[nsets];
TGraphErrors *gVnEP[nsets];
TGraphErrors *gVnSP[nsets];

TGraphErrors *gRnonuni[nsets];
TGraphErrors *gRtrueNonuni[nsets];
TGraphErrors *gVnNonuni[nsets];
TGraphErrors *gVnTrueNonuni[nsets];
TGraphErrors *gVnEPnonuni[nsets];
TGraphErrors *gVnSPnonuni[nsets];

TGraphErrors *gPtBin[nsets];

int mMarker[nsets] = {24};

double mSize = 1.0;

TFile *fIn[nsets];

double VnDist(double *x, double *p);
void hset(TH1& hid, TString xtit="", TString ytit="",
		double titoffx = 1.1, double titoffy = 1.1,
		double titsizex = 0.06, double titsizey = 0.06,
		double labeloffx = 0.01, double labeloffy = 0.001,
		double labelsizex = 0.05, double labelsizey = 0.05,
		int divx = 505, int divy=505);

void PlotVn() {

    TH2F *hfr;

    int i, j;
    for (i=0; i<nsets; i++) {
        fIn[i] = TFile::Open(fileNames[i], "read");
        if(fIn[i]==0) ErrorExit(Form("Cannot open file: %s",fileNames[i].Data()));

        hInputFlow[i] = (TH1D*) fIn[i]->Get("hInputFlow");

        gR[i] = (TGraphErrors*) fIn[i]->Get("gR");
        gRtrue[i] = (TGraphErrors*) fIn[i]->Get("gRtrue");
        gRnonuni[i] = (TGraphErrors*) fIn[i]->Get("gRnonuni");
        gRtrueNonuni[i] = (TGraphErrors*) fIn[i]->Get("gRtrueNonuni");

        gVn[i] = (TGraphErrors*) fIn[i]->Get("gVn");
        gVnTrue[i] = (TGraphErrors*) fIn[i]->Get("gVnTrue");
        gVnEP[i] = (TGraphErrors*) fIn[i]->Get("gVnEP");
        gVnSP[i] = (TGraphErrors*) fIn[i]->Get("gVnSP");

        gVnNonuni[i] = (TGraphErrors*) fIn[i]->Get("gVnNonuni");
        gVnTrueNonuni[i] = (TGraphErrors*) fIn[i]->Get("gVnTrueNonuni");
        gVnEPnonuni[i] = (TGraphErrors*) fIn[i]->Get("gVnEPnonuni");
        gVnSPnonuni[i] = (TGraphErrors*) fIn[i]->Get("gVnSPnonuni");

        gPtBin[i] = (TGraphErrors*) fIn[i]->Get("gPtBin");

        gR[i]->SetMarkerStyle(mMarker[i]+1);
        gR[i]->SetMarkerColor(i+1);
        gR[i]->SetMarkerSize(mSize);

        gRtrue[i]->SetMarkerStyle(mMarker[i]+1);
        gRtrue[i]->SetMarkerColor(i+2);
        gRtrue[i]->SetMarkerSize(mSize);

        gRnonuni[i]->SetMarkerStyle(mMarker[i]+1);
        gRnonuni[i]->SetMarkerColor(i+3);
        gRnonuni[i]->SetMarkerSize(mSize);

        gRtrueNonuni[i]->SetMarkerStyle(mMarker[i]+1);
        gRtrueNonuni[i]->SetMarkerColor(i+4);
        gRtrueNonuni[i]->SetMarkerSize(mSize);

        gVn[i]->SetMarkerStyle(mMarker[i]);
        gVn[i]->SetMarkerColor(i+1);
        gVn[i]->SetMarkerSize(mSize);

        gVnTrue[i]->SetMarkerStyle(mMarker[i]);
        gVnTrue[i]->SetMarkerColor(i+2);
        gVnTrue[i]->SetMarkerSize(mSize);

        gVnEP[i]->SetMarkerStyle(mMarker[i]);
        gVnEP[i]->SetMarkerColor(i+3);
        gVnEP[i]->SetMarkerSize(mSize);

        gVnSP[i]->SetMarkerStyle(mMarker[i]);
        gVnSP[i]->SetMarkerColor(i+4);
        gVnSP[i]->SetMarkerSize(mSize);

        gVnNonuni[i]->SetMarkerStyle(mMarker[i]);
        gVnNonuni[i]->SetMarkerColor(i+1);
        gVnNonuni[i]->SetMarkerSize(mSize);

        gVnTrueNonuni[i]->SetMarkerStyle(mMarker[i]);
        gVnTrueNonuni[i]->SetMarkerColor(i+2);
        gVnTrueNonuni[i]->SetMarkerSize(mSize);

        gVnEPnonuni[i]->SetMarkerStyle(mMarker[i]);
        gVnEPnonuni[i]->SetMarkerColor(i+3);
        gVnEPnonuni[i]->SetMarkerSize(mSize);

        gVnSPnonuni[i]->SetMarkerStyle(mMarker[i]);
        gVnSPnonuni[i]->SetMarkerColor(i+4);
        gVnSPnonuni[i]->SetMarkerSize(mSize);

        gPtBin[i]->SetMarkerStyle(mMarker[i]);
        gPtBin[i]->SetMarkerColor(i+2);
        gPtBin[i]->SetMarkerSize(mSize);

    }

    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1", "vn values with different methods");
    c1->cd();

    TLegend *leg1 = new TLegend(0.50,0.6,0.85,0.85,"","brNDC");
    leg1->SetTextSize(0.037);leg1->SetBorderSize(0);

    hfr = new TH2F(Form("hfr1%d",0)," ", 1, 0.5, 5.5, 1, -0.021, 0.18);
    hset( *hfr, "n", "v_{n}",1.0,1.0, 0.04,0.04, 0.01,0.01, 0.03,0.03, 510,510);
    hfr->Draw();

    for (i=0; i<nsets; i++) {
        if (i==0) {
            gVnTrue[i]->Draw("SAME P");
            gVnTrue[i]->SetTitle("Vn values with different methods - uniform phi; Vn; n");
            leg1->AddEntry(gVnTrue[i], "v_{n}, true RP", "p");
        }

        gVn[i]->Draw("SAME P");
        leg1->AddEntry(gVn[i], "v_{n}, trad. EP", "p");

        gVnEP[i]->Draw("SAME P");
        leg1->AddEntry(gVnEP[i], "v_{n}{EP}", "p");

        gVnSP[i]->Draw("SAME P");
        leg1->AddEntry(gVnSP[i], "v_{n}{SP}", "p");

        hInputFlow[i]->Draw("SAME HIST");
        leg1->AddEntry(hInputFlow[i], "Input", "l");
    }

    leg1->Draw("SAME");
    c1->SaveAs("figures/vn.pdf");

    TCanvas *c2 = new TCanvas("c2", "vn values with different methods - nonuniform phi");
    c2->cd();

    TLegend *leg2 = new TLegend(0.50,0.6,0.85,0.85,"","brNDC");
    leg2->SetTextSize(0.037);leg2->SetBorderSize(0);

    hfr = new TH2F(Form("hfr2%d",0)," ", 1, 0.5, 5.5, 1, -0.021, 0.18);
    hset( *hfr, "n", "v_{n}",1.0,1.0, 0.04,0.04, 0.01,0.01, 0.03,0.03, 510,510);
    hfr->Draw();

    for (i=0; i<nsets; i++) {
        if (i==0) {
            gVnTrueNonuni[i]->Draw("SAME P");
            gVnTrueNonuni[i]->SetTitle("Vn values with different methods; Vn; n");
            leg2->AddEntry(gVnTrueNonuni[i], "v_{n}, true RP", "p");
        }

        gVnNonuni[i]->Draw("SAME P");
        leg2->AddEntry(gVnNonuni[i], "v_{n}, trad. EP", "p");

        gVnEPnonuni[i]->Draw("SAME P");
        leg2->AddEntry(gVnEPnonuni[i], "v_{n}{EP}", "p");

        gVnSPnonuni[i]->Draw("SAME P");
        leg2->AddEntry(gVnSPnonuni[i], "v_{n}{SP}", "p");

        hInputFlow[i]->Draw("SAME HIST");
        leg2->AddEntry(hInputFlow[i], "Input", "l");
    }

    leg2->Draw("SAME");
    c2->SaveAs("figures/vn-nonuni.pdf");

    TCanvas *c3 = new TCanvas("c3", "R values with different methods");
    c3->cd();

    hfr = new TH2F(Form("hfr%d",1)," ", 1, 0.5, 5.5, 1, -0.021, 1.04);
    hset( *hfr, "n", "R_{n}",1.0,1.0, 0.04,0.04, 0.01,0.01, 0.03,0.03, 510,510);
    hfr->Draw();

    TLegend *leg3 = new TLegend(0.25,0.15,0.55,0.30,"","brNDC");
    leg3->SetTextSize(0.037);leg3->SetBorderSize(0);

    for (i=0; i<nsets; i++) {
        if (i==0) {
            gRtrue[i]->Draw("SAME P");
            gRtrue[i]->SetTitle("Resolution parameter; R; n");
            leg3->AddEntry(gRtrue[i], "R true", "p");
        }

        gRtrueNonuni[i]->Draw("SAME P");
        leg3->AddEntry(gRtrueNonuni[i], "R true, nonuni", "p");

        gR[i]->Draw("SAME P");
        leg3->AddEntry(gR[i], "R sub event method", "p");

        gRnonuni[i]->Draw("SAME P");
        leg3->AddEntry(gRnonuni[i], "R sub event method, nonuni", "p");
    }

    leg3->Draw("SAME");
    c3->SaveAs("figures/Rn.pdf");

    TCanvas *c4 = new TCanvas("c4", "v2(pT)");
    c4->cd();

    TF1 *fVnDist = new TF1("fVnDist", VnDist, 0.0, 10.0, 3);
    fVnDist->SetParameters(2.0, 1.0, 0.15);

    for (i=0; i<nsets; i++) {
        if (i==0) {
            gPtBin[i]->Draw("AP");
            gPtBin[i]->SetTitle("v2 as function of pT; v2; pT");
        }

        fVnDist->Draw("SAME");
    }
    c4->SaveAs("figures/v2Pt.pdf");

    TH1D *hPhi = (TH1D*)fIn[0]->Get("hPhi");
    TH1D *hPhiNonuni = (TH1D*)fIn[0]->Get("hPhiNonuni");

    TCanvas *c5 = new TCanvas("c5", "c5");
    c5->Divide(2,1);

    c5->cd(1);
    hPhi->Draw("HIST");

    c5->cd(2);
    hPhiNonuni->Draw("HIST");

    c5->SaveAs("figures/phi.pdf");
}

double VnDist(double *x, double *p) {
    double pt = x[0];
    double alpha = p[0];
    double beta = p[1];
    double vnMax = p[2];
    double C = vnMax/(TMath::Power(alpha/beta, alpha)*TMath::Exp(-alpha));
    return C*TMath::Power(pt, alpha)*TMath::Exp(-beta*pt);
}

void hset(TH1& hid, TString xtit="", TString ytit="",
		double titoffx = 1.1, double titoffy = 1.1,
		double titsizex = 0.06, double titsizey = 0.06,
		double labeloffx = 0.01, double labeloffy = 0.001,
		double labelsizex = 0.05, double labelsizey = 0.05,
		int divx = 505, int divy=505)
{
	hid.GetXaxis()->CenterTitle(1);
	hid.GetYaxis()->CenterTitle(1);

	hid.GetXaxis()->SetTitleOffset(titoffx);
	hid.GetYaxis()->SetTitleOffset(titoffy);

	hid.GetXaxis()->SetTitleSize(titsizex);
	hid.GetYaxis()->SetTitleSize(titsizey);

	hid.GetXaxis()->SetLabelOffset(labeloffx);
	hid.GetYaxis()->SetLabelOffset(labeloffy);

	hid.GetXaxis()->SetLabelSize(labelsizex);
	hid.GetYaxis()->SetLabelSize(labelsizey);

	hid.GetXaxis()->SetNdivisions(divx);
	hid.GetYaxis()->SetNdivisions(divy);

	hid.GetXaxis()->SetTitle(xtit);
	hid.GetYaxis()->SetTitle(ytit);

	hid.GetXaxis()->SetLabelFont(42);
	hid.GetYaxis()->SetLabelFont(42);
	hid.GetXaxis()->SetTitleFont(42);
	hid.GetYaxis()->SetTitleFont(42);
}
