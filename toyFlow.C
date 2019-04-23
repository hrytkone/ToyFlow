#include <iostream>
#include <stdlib.h>

#include "JHistos.h"
#include "JToyMCTrack.h"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include "TComplex.h"
#include "TCanvas.h"

#include "TStopwatch.h"

using namespace std;

double PtDist(double *x, double *p);
double PhiDist(double *x, double *p);
double VnDist(double *x, double *p);

void GetEvent(TClonesArray *listUni, TClonesArray *listUniA, TClonesArray *listUniB, TClonesArray *listNonuni, TClonesArray *listNonuniA, TClonesArray *listNonuniB, JHistos *histos, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, TF1 *fVnDist, double *vn, double *Psi, double alpha, double beta, int nMult, double percentage, double phiMin, double phiMax, bool bUsePtDependence);
void AnalyzeEvent(TClonesArray *listFull, TClonesArray *listSubA, TClonesArray *listSubB, JHistos *histos, double *Psi, bool bUseWeight, bool bDoCorrections, double **corrections);

double AcceptanceFunc(double *x, double *p);
double AcceptanceFuncTimesSin(double *x, double *p);
double AcceptanceFuncTimesCos(double *x, double *p);
double *GetCorrectionParam(double phiMin, double phiMax, double percentage, double n, TF1 *fA, TF1 *fTimesSin, TF1 *fTimesCos, TF1 *fTimesSin2, TF1 *fTimesCos2);

double ShiftCorrection(double x, double correction);
double TwistCorrection(double x, double y, double lambda1, double lambda2);
double RescalingCorrection(double x, double a);
void DoCorrections(TComplex &Q, double cm, double sm, double lambdaMinus, double lambdaPlus, double aMinus, double aPlus);

double GetEventPlane(TComplex Q, int n);
double GetVnObs(TComplex Q, double phi, int n);
double GetScalarProduct(TComplex Qa, TComplex Qb);

double CheckIfZero(double x, double thres);

int main(int argc, char **argv) {

    TString outFileName = argc > 1 ? argv[1]:"toyFlow.root";
    int nEvents = argc > 2 ? atol(argv[2]):1000;
    const double dNdeta = argc > 3 ? atof(argv[3]):1500.0;
    const double etaRange = 0.8;
    const int nMult = int(2. * etaRange * dNdeta);
    cout << "dNdeta=" << dNdeta << "     total multiplicity=" << nMult << endl;

    const int nCoef = 5;

    double alpha = 2.0;
    double beta = 1.0;

    const double scale = 1.0;
    double vn[nCoef] = {scale*0.0, scale*0.15, scale*0.08, scale*0.03, scale*0.01};

    bool bRandomPsi = true;
    bool bUsePtDependence = false;
    bool bUseWeight = false;

    int i, j; // indices for loops
    double pi = TMath::Pi();

    TStopwatch timer;
    timer.Start();

    TFile *fOut = TFile::Open(outFileName, "RECREATE");

    TRandom3 *rand = new TRandom3(0);
    rand->SetSeed(0);

    JHistos *histos = new JHistos();

    double Psi[nCoef] = {0};

    double Tdec = 0.12;
    double vr = 0.6;
    double Teff = Tdec * TMath::Sqrt((1.+vr)/(1.-vr));
    TF1 *fPtDist = new TF1("fPtDist", PtDist, 0.0, 10.0, 1);
    fPtDist->SetParameter(0, 1./Teff);

    TF1 *fPhiDist = new TF1("fPhiDist", PhiDist, -pi, pi, 10);
    fPhiDist->SetParameters(vn[0], vn[1], vn[2], vn[3], vn[4], Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

    TF1 *fVnDist = new TF1("fVnDist", VnDist, 0.0, 10.0, 3);

    TClonesArray *listUniFull = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *listUniA = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *listUniB = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *listNonuniFull = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *listNonuniA = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *listNonuniB = new TClonesArray("JToyMCTrack", nMult+1);

    // CORRECTIONS
    double phiMin = 1.0;
    double phiMax = 1.5;
    double percentage = 30.0;

    TF1 *fA = new TF1("fA", AcceptanceFunc, -pi, pi, 4);
    TF1 *fTimesSin = new TF1("fTimesSin", AcceptanceFuncTimesSin, -pi, pi, 5);
    TF1 *fTimesCos = new TF1("fTimesCos", AcceptanceFuncTimesCos, -pi, pi, 5);
    TF1 *fTimesSin2 = new TF1("fTimesSin2", AcceptanceFuncTimesSin, -pi, pi, 5);
    TF1 *fTimesCos2 = new TF1("fTimesCos2", AcceptanceFuncTimesCos, -pi, pi, 5);

    cout << "CORRECTIONS:\n";
    cout << "vn: {cm, sm, cm2, sm2, a-, a+, lambda-, lambda+}\n";

    double *corrections[5];
    for (i=0; i<5; i++) {
        corrections[i] = GetCorrectionParam(phiMin, phiMax, percentage, i+1, fA, fTimesSin, fTimesCos, fTimesSin2, fTimesCos2);
    }

    // Save correction parameters
    TH2D *hCorrectionParameters = new TH2D("hCorrectionParameters", "hCorrectionParameters", 5, 0.5, 5.5, 8, 0.5, 8.5);
    for (i=0; i<5; i++) {
        for (j=0; j<8; j++) {
            hCorrectionParameters->Fill(i, j, corrections[i][j]);
        }
    }

    // Save input numbers
    TH1D *hInputNumbers = new TH1D("hInputNumbers","hInputNumbers",16, 0.5, 16.5);
    hInputNumbers->Fill(1, double(nEvents));
    hInputNumbers->Fill(2, double(dNdeta));
    hInputNumbers->Fill(3, etaRange);
    hInputNumbers->Fill(4, double(nMult));
    hInputNumbers->Fill(5, vn[0]);
    hInputNumbers->Fill(6, vn[1]);
    hInputNumbers->Fill(7, vn[2]);
    hInputNumbers->Fill(8, vn[3]);
    hInputNumbers->Fill(9, vn[4]);
    hInputNumbers->Fill(10, Tdec);
    hInputNumbers->Fill(11, vr);
    hInputNumbers->Fill(12, Teff);
    hInputNumbers->Fill(13, 1./Teff);
    hInputNumbers->Fill(14, phiMin);
    hInputNumbers->Fill(15, phiMax);
    hInputNumbers->Fill(16, percentage);

    int nOutput = nEvents/20;
    for (i=0; i<nEvents; i++) {
        if (i % nOutput == 0)
            cout << 100*i/nEvents << " % finished" << endl;

        // Good to use argument "C" in case of pointers
        listUniFull->Clear("C");
        listUniA->Clear("C");
        listUniB->Clear("C");
        listNonuniFull->Clear("C");
        listNonuniA->Clear("C");
        listNonuniB->Clear("C");

        if (bRandomPsi) {
            for (j=0; j<5; j++) {
                Psi[j] = rand->Uniform(0, 2*pi);
            }
        } else {
            double psiTemp = rand->Uniform(0, 2*pi);
            for (j=0; j<5;j++) {
                Psi[j] = psiTemp;
            }
        }

        fPhiDist->SetParameters(vn[0], vn[1], vn[2], vn[3], vn[4],
            Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

        GetEvent(listUniFull, listUniA, listUniB, listNonuniFull, listNonuniA, listNonuniB, histos, rand, fPtDist, fPhiDist, fVnDist, vn, Psi, alpha, beta, nMult, percentage, phiMin, phiMax, bUsePtDependence);
        AnalyzeEvent(listUniFull, listUniA, listUniB, histos, Psi, bUseWeight, false, corrections);
        AnalyzeEvent(listNonuniFull, listNonuniA, listNonuniB, histos, Psi, bUseWeight, true, corrections);

    }

    fOut->Write();
    timer.Print();
    return 0;
}

//======END OF MAIN PROGRAM======
void GetEvent(TClonesArray *listUni, TClonesArray *listUniA, TClonesArray *listUniB, TClonesArray *listNonuni, TClonesArray *listNonuniA, TClonesArray *listNonuniB, JHistos *histos, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, TF1 *fVnDist, double *vn, double *Psi, double alpha, double beta, int nMult, double percentage, double phiMin, double phiMax, bool bUsePtDependence) {

    double pT, phi, eta, E;
    double px, py, pz;
    int nNonuniFull = 0;
    int nUniA = 0;
    int nUniB = 0;
    int nNonuniA = 0;
    int nNonuniB = 0;

    JToyMCTrack track;
    TLorentzVector lVec;

    double randNum;
    double vnTemp[5];
    int i, j, m;
    for (i = 0; i < nMult; i++) {

        pT = fPt->GetRandom();
        if (bUsePtDependence) {
            for (j=0; j<5; j++) {
                fVnDist->SetParameters(alpha, beta, vn[j]);
                vnTemp[j] = fVnDist->Eval(pT);
            }
            fPhi->SetParameters(vnTemp[0], vnTemp[1], vnTemp[2], vnTemp[3], vnTemp[4],
                Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);
        }

        phi = fPhi->GetRandom();
        eta = rand->Uniform(-0.8, 0.8);

        histos->hPt->Fill(pT);
        histos->hPhi->Fill(phi);

        px = pT*TMath::Cos(phi);
        py = pT*TMath::Sin(phi);
        pz = pT*TMath::SinH(eta);
        E = TMath::Sqrt(pT*pT + pz*pz);
        lVec.SetPxPyPzE(px, py, pz, E);
        track.SetTrack(lVec);

        new((*listUni)[i]) JToyMCTrack(track);
        if (i%2 == 0) {
            new((*listUniA)[nUniA]) JToyMCTrack(track);
            nUniA++;
        } else {
            new((*listUniB)[nUniB]) JToyMCTrack(track);
            nUniB++;
        }

        if (phi < phiMin || phi > phiMax) {

            histos->hPhiNonuni->Fill(phi);

            new((*listNonuni)[nNonuniFull]) JToyMCTrack(track);
            nNonuniFull++;
            if ( nNonuniA<nNonuniB ) {
                new((*listNonuniA)[nNonuniA]) JToyMCTrack(track);
                nNonuniA++;
            } else {
                new((*listNonuniB)[nNonuniB]) JToyMCTrack(track);
                nNonuniB++;
            }
        } else {
            randNum = rand->Rndm();
            if ( randNum > percentage ) {

                histos->hPhiNonuni->Fill(phi);

                new((*listNonuni)[nNonuniFull]) JToyMCTrack(track);
                nNonuniFull++;
                if ( nNonuniA<nNonuniB ) {
                    new((*listNonuniA)[nNonuniA]) JToyMCTrack(track);
                    nNonuniA++;
                } else {
                    new((*listNonuniB)[nNonuniB]) JToyMCTrack(track);
                    nNonuniB++;
                }

            }
        }

        histos->hMultiplicity->Fill(nMult);
        histos->hMultiplicityNonuni->Fill(nNonuniFull);
    }
}

void AnalyzeEvent(TClonesArray *listFull, TClonesArray *listSubA, TClonesArray *listSubB, JHistos *histos, double *Psi, bool bUseWeight, bool bDoCorrections, double **corrections) {

    int nMult = listFull->GetEntriesFast();
    int nMultPOI = listSubA->GetEntriesFast();
    int nMultRef = listSubB->GetEntriesFast();

    int i, j, k, n;
    double w = 1.0;
    double phi, pt;
    double Rtrue;

    TComplex Qfull = TComplex(0, 0);
    TComplex QPOI = TComplex(0, 0);
    TComplex Qref = TComplex(0, 0);
    TComplex QrefA = TComplex(0, 0);
    TComplex QrefB = TComplex(0, 0);
    TComplex unitVec = TComplex(0, 0);
    TComplex autocorr = TComplex(0, 0);

    double cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus;
    double EventPlaneA, EventPlaneB, Rsub, vobs;

    double normFull, normPOI, normRef, normA, normB;

    double QnQnRef, QnAQnB;

    int nPtBins = 9;
    double point = 0.0;
    double step = 0.0;
    double binWidths[nPtBins+1];
    for (i=0; i<nPtBins+1; i++) {
        binWidths[i] = point;
        step += 0.2;
        point += step;
    }
    vector<vector<TComplex>> pTBinsQ;
    pTBinsQ.resize(nPtBins);

    JToyMCTrack *track;

    for (i = 0; i < 5; i++) {

        n = i+1;

        cm = corrections[i][0];
        sm = corrections[i][1];
        aMinus = corrections[i][4];
        aPlus = corrections[i][5];
        lambdaMinus = corrections[i][6];
        lambdaPlus = corrections[i][7];

        Qfull = TComplex(0, 0);
        QPOI = TComplex(0, 0);
        Qref = TComplex(0, 0);
        QrefA = TComplex(0, 0);
        QrefB = TComplex(0, 0);

        vobs = 0.0;

        normFull = 0.0; normPOI = 0.0; normRef = 0.0; normA = 0.0; normB = 0.0;

        for (j=0; j<nMult; j++) {

            track = (JToyMCTrack*)listFull->At(j);
            phi = track->GetPhi();

            if (bUseWeight) w = track->GetPt();
            normFull += w*w;

            unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
            if (bDoCorrections) DoCorrections(unitVec, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
            Qfull += unitVec;

            if ( j<nMultPOI ) {
                track = (JToyMCTrack*)listSubA->At(j);
                phi = track->GetPhi();

                if (bUseWeight) w = track->GetPt();
                normPOI += w*w;

                unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
                if (bDoCorrections) DoCorrections(unitVec, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
                QPOI += unitVec;
            }

            if ( j<nMultRef ) {
                track = (JToyMCTrack*)listSubB->At(j);
                phi = track->GetPhi();

                if (bUseWeight) w = track->GetPt();
                normRef += w*w;

                unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
                if (bDoCorrections) DoCorrections(unitVec, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
                Qref += unitVec;

                if (j<nMultRef/2.0) {
                    normA += w*w;
                    QrefA += unitVec;
                } else {
                    normB += w*w;
                    QrefB += unitVec;
                }
            }
        }

        for (j = 0; j < nMult; j++) {

            track = (JToyMCTrack*)listFull->At(j);
            phi = track->GetPhi();

            autocorr = TComplex(TMath::Cos(n*phi), TMath::Sin(n*phi));

            Qfull -= autocorr;
            vobs += GetVnObs(Qfull, phi, n);
            Qfull += autocorr;
        }

        vobs /= nMult;

        // Resolution parameter calculations
        Rtrue = TMath::Cos(n*(GetEventPlane(Qfull, n) - Psi[i]));

        EventPlaneA = GetEventPlane(QrefA, n);
        EventPlaneB = GetEventPlane(QrefB, n);
        Rsub = TMath::Cos(n*(EventPlaneA - EventPlaneB));

        // ALTERNATIVE EVENT PLANE METHOD
        normPOI = TMath::Sqrt(normPOI);
        normRef = TMath::Sqrt(normRef);
        normA = TMath::Sqrt(normA);
        normB = TMath::Sqrt(normB);

        QPOI /= normPOI; Qref /= normRef; QrefA /= normA; QrefB /= normB;

        QnQnRef = GetScalarProduct(QPOI, Qref);
        QnQnRef /= TComplex::Abs(Qref);

        QnAQnB = GetScalarProduct(QrefA, QrefB);
        QnAQnB /= TComplex::Abs(QrefA);
        QnAQnB /= TComplex::Abs(QrefB);

        if (bDoCorrections) {
            histos->hVnObsCorrected[i]->Fill(vobs);
            histos->hRtrueCorrected[i]->Fill(Rtrue);
            histos->hRsubCorrected[i]->Fill(Rsub);
            histos->hQnQnAcorrected[i]->Fill(QnQnRef);
            histos->hQnAQnBcorrected[i]->Fill(QnAQnB);
        } else {
            histos->hVnObs[i]->Fill(vobs);
            histos->hRtrue[i]->Fill(Rtrue);
            histos->hRsub[i]->Fill(Rsub);
            histos->hQnQnA[i]->Fill(QnQnRef);
            histos->hQnAQnB[i]->Fill(QnAQnB);
        }

        // Divide into pT-bins
        if (n==2 && !bDoCorrections) {
            double weight = 0.0;
            double norms[nPtBins];
            for (j=0; j<nPtBins; j++) norms[j] = 0;

            for (j=0; j < nMultPOI; j++) {

                track = (JToyMCTrack*)listFull->At(j);
                phi = track->GetPhi();
                pt = track->GetPt();

                if (bUseWeight) w = pt;

                unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));

                for (k=0; k<nPtBins; k++) {
                    if ((binWidths[k] <= pt) && (binWidths[k+1] > pt)) {
                        pTBinsQ[k].push_back(unitVec);
                        norms[k] += w*w;
                    }
                }
            }

            for (j=0; j<nPtBins; j++) {
                QPOI = TComplex(0, 0);
                for (k=0; k<pTBinsQ[j].size(); k++) {
                    QPOI += pTBinsQ[j][k];
                }
                weight = TMath::Sqrt(norms[j]);
                QPOI /= weight;
                QnQnRef = GetScalarProduct(QPOI, Qref);
                QnQnRef /= TComplex::Abs(Qref);
                histos->hPtBin[j]->Fill(QnQnRef);
                histos->hSqrtSumWeightsPtBins[j]->Fill(weight);
                pTBinsQ[j].clear();
            }

        }
    }

    if (bDoCorrections) {
        histos->hSqrtSumWeightsNonuni->Fill(normPOI);
        histos->hSqrtSumWeightsANonuni->Fill(normA);
        histos->hSqrtSumWeightsBNonuni->Fill(normB);
    } else {
        histos->hSqrtSumWeights->Fill(normPOI);
        histos->hSqrtSumWeightsA->Fill(normA);
        histos->hSqrtSumWeightsB->Fill(normB);
    }
}

double PtDist(double *x, double *p) {
    return TMath::Exp(-p[0]*x[0]);
}

double PhiDist(double *x, double *p) {
    double phi = x[0];
    double vn[5] = {p[0], p[1], p[2], p[3], p[4]};
    double psi[5] = {p[5], p[6], p[7], p[8], p[9]};
    return 1.0 + 2.0*vn[0]*TMath::Cos(phi - psi[0]) + 2.0*vn[1]*TMath::Cos(2.*(phi - psi[1]))
        + 2.0*vn[2]*TMath::Cos(3.*(phi - psi[2])) + 2.0*vn[3]*TMath::Cos(4.*(phi - psi[3]))
        + 2.0*vn[4]*TMath::Cos(5.*(phi - psi[4]));
}

double VnDist(double *x, double *p) {
    double pt = x[0];
    double alpha = p[0];
    double beta = p[1];
    double vnMax = p[2];

    double C = vnMax/(TMath::Power(alpha/beta, alpha)*TMath::Exp(-alpha));
    return C*TMath::Power(pt, alpha)*TMath::Exp(-beta*pt);
}

double AcceptanceFunc(double *x, double *p) {
    double pi = TMath::Pi();
    double phi = x[0];
    double phiMin = p[0];
    double phiMax = p[1];
    double missedParticles = p[2];
    double norm = p[3];
    return ((phi < phiMin) || (phi > phiMax)) ? norm/(2*pi) : norm*(1-missedParticles)/(2*pi);
}

double AcceptanceFuncTimesSin(double *x, double *p) {
    double pi = TMath::Pi();
    double phi = x[0];
    double phiMin = p[0];
    double phiMax = p[1];
    double missedParticles = p[2];
    double n = p[3];
    double norm = p[4];
    return ((phi < phiMin) || (phi > phiMax)) ? norm*TMath::Sin(n*phi)/(2*pi) : norm*(1-missedParticles)*TMath::Sin(n*phi)/(2*pi);
}

double AcceptanceFuncTimesCos(double *x, double *p) {
    double pi = TMath::Pi();
    double phi = x[0];
    double phiMin = p[0];
    double phiMax = p[1];
    double missedParticles = p[2];
    double n = p[3];
    double norm = p[4];
    return ((phi < phiMin) || (phi > phiMax)) ? norm*TMath::Cos(n*phi)/(2*pi) : norm*(1-missedParticles)*TMath::Cos(n*phi)/(2*pi);
}

double *GetCorrectionParam(double phiMin, double phiMax, double percentage, double n, TF1 *fA, TF1 *fTimesSin, TF1 *fTimesCos, TF1 *fTimesSin2, TF1 *fTimesCos2) {
    double pi = TMath::Pi();
    fA->SetParameters(phiMin, phiMax, percentage, 1);
    double area = fA->Integral(-pi, pi);

    fTimesSin->SetParameters(phiMin, phiMax, percentage, n, 1.0/area);
    fTimesCos->SetParameters(phiMin, phiMax, percentage, n, 1.0/area);
    fTimesSin2->SetParameters(phiMin, phiMax, percentage, 2*n, 1.0/area);
    fTimesCos2->SetParameters(phiMin, phiMax, percentage, 2*n, 1.0/area);

    double thres = 0.0001;
    double *corrections = new double[8];
    corrections[0] = CheckIfZero(fTimesCos->Integral(-pi, pi), thres); //cm
    corrections[1] = CheckIfZero(fTimesSin->Integral(-pi, pi), thres); //sm
    corrections[2] = CheckIfZero(fTimesCos2->Integral(-pi, pi), thres); //cm2
    corrections[3] = CheckIfZero(fTimesSin2->Integral(-pi, pi), thres); //sm2
    corrections[4] = CheckIfZero(1.0 - corrections[2], thres); //a-
    corrections[5] = CheckIfZero(1.0 + corrections[2], thres); //a+
    corrections[6] = CheckIfZero(corrections[3]/corrections[4], thres); //lambda-
    corrections[7] = CheckIfZero(corrections[3]/corrections[5], thres); //lambda+

    cout << "v" << n << ": {";
    for (int i=0; i<7; i++) {
        cout << corrections[i] << ", ";
    }
    cout << corrections[7] << "}\n";

    return corrections;
}

double CheckIfZero(double x, double thres) {
    if (abs(x) < thres) {
        return 0.0;
    } else {
        return x;
    }
}

double ShiftCorrection(double x, double correction) {
    return x - correction;
}

double TwistCorrection(double x, double y, double lambda1, double lambda2) {
    return (x - lambda1*y)/(1-lambda1*lambda2);
}

double RescalingCorrection(double x, double a) {
    return x/a;
}

void DoCorrections(TComplex &Q, double cm, double sm, double lambdaMinus, double lambdaPlus, double aMinus, double aPlus) {
    Q = TComplex(ShiftCorrection(Q.Re(), cm), ShiftCorrection(Q.Im(), sm));
    Q = TComplex(TwistCorrection(Q.Re(), Q.Im(), lambdaMinus, lambdaPlus), TwistCorrection(Q.Im(), Q.Re(), lambdaPlus, lambdaMinus));
    Q = TComplex(RescalingCorrection(Q.Re(), aPlus), RescalingCorrection(Q.Im(), aMinus));
}

double GetEventPlane(TComplex Q, int n) {
    double Qx = Q.Re();
    double Qy = Q.Im();
    return TMath::ATan2(Qy, Qx)/n;
}

double GetVnObs(TComplex Q, double phi, int n) {
    double EP = GetEventPlane(Q, n);
    return TMath::Cos(n*(phi - EP));
}

double GetScalarProduct(TComplex Qa, TComplex Qb) {
    return Qa.Re()*Qb.Re() + Qa.Im()*Qb.Im();
}
