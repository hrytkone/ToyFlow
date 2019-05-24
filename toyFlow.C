#include <iostream>
#include <stdlib.h>

// OWN
#include "src/JHistos.h"
#include "src/JEventLists.h"
#include "src/JToyMCTrack.h"
#include "src/JInputs.h"
#include "src/JConst.h"

// ROOT
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include "TComplex.h"
#include "TCanvas.h"

// OTHER
#include "TStopwatch.h"

using namespace std;

double PtDist(double *x, double *p);
double PhiDist(double *x, double *p);
double VnDist(double *x, double *p);

void GetEvent(JHistos *histos, JEventLists *lists, JInputs *inputs, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, TF1 *fVnDist, double *vn, double *Psi, double percentage, double phiMin, double phiMax, bool bUsePtDependence, bool bUseGranularity, double startAngle);
void AnalyzeEvent(JEventLists *lists, JHistos *histos, double *Psi, bool bUseWeight, bool bDoCorrections, double **corrections);

double AcceptanceFunc(double *x, double *p);
double AcceptanceFuncTimesSin(double *x, double *p);
double AcceptanceFuncTimesCos(double *x, double *p);
double *GetCorrectionParam(double phiMin, double phiMax, double percentage, double n, TF1 *fA, TF1 *fTimesSin, TF1 *fTimesCos, TF1 *fTimesSin2, TF1 *fTimesCos2);

double ShiftCorrection(double x, double correction);
double TwistCorrection(double x, double y, double lambda1, double lambda2);
double RescalingCorrection(double x, double a);
void DoCorrections(TComplex &Qvec, double cm, double sm, double lambdaMinus, double lambdaPlus, double aMinus, double aPlus);

void CalculateQvector(JToyMCTrack *track, TComplex unitVec, TComplex &Qvec, double &norm, bool bUseWeight, bool bDoCorrections, int n, double w, double cm, double sm, double lambdaMinus, double lambdaPlus, double aMinus, double aPlus);
double GetEventPlane(TComplex Qvec, int n);
double GetVnObs(TComplex Qvec, double phi, int n);

double CheckIfZero(double x, double thres);
double CheckPhi(double phi, double startAngle);
double CheckEta(double eta);

int main(int argc, char **argv) {

    TString outFileName = argc > 1 ? argv[1]:"toyFlow.root";
    int nEvents = argc > 2 ? atol(argv[2]):1000;
    bool bUsePtDependence = argc > 3 ? atol(argv[3]):0;
    bool bUseGranularity = argc > 4 ? atol(argv[4]):0;
    int iSeed = argc > 5 ? atol(argv[5]):0;

    bool bUseWeight = false;
    bool bRandomPsi = true;

    const double scale = 1.0;
    double vn[nCoef] = {scale*0.0, scale*0.15, scale*0.08, scale*0.03, scale*0.01};

    int i, j; // indices for loops

    TStopwatch timer;
    timer.Start();

    TFile *fOut = TFile::Open(outFileName, "RECREATE");

    TRandom3 *rand = new TRandom3(iSeed);
    //rand->SetSeed(0);

    JHistos *histos = new JHistos();
    JEventLists *lists = new JEventLists();
    JInputs *inputs = new JInputs();
    inputs->Load();

    double Psi[nCoef] = {0};

    double Teff = Tdec * TMath::Sqrt((1.+vr)/(1.-vr));
    TF1 *fPtDist = new TF1("fPtDist", PtDist, 0.0, 10.0, 1);
    fPtDist->SetParameter(0, 1./Teff);

    TF1 *fPhiDist = new TF1("fPhiDist", PhiDist, -PI, PI, 10);
    fPhiDist->SetParameters(vn[0], vn[1], vn[2], vn[3], vn[4], Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

    TF1 *fVnDist = new TF1("fVnDist", VnDist, 0.0, 10.0, 3);

    // CORRECTIONS
    double phiMin = 1.12;
    double phiMax = 1.54;
    double percentage = 0.5; // REMEMBER TO GIVE PERCENTAGE IN RANGE 0.0-1.0!

    TF1 *fA = new TF1("fA", AcceptanceFunc, -PI, PI, 4);
    TF1 *fTimesSin = new TF1("fTimesSin", AcceptanceFuncTimesSin, -PI, PI, 5);
    TF1 *fTimesCos = new TF1("fTimesCos", AcceptanceFuncTimesCos, -PI, PI, 5);
    TF1 *fTimesSin2 = new TF1("fTimesSin2", AcceptanceFuncTimesSin, -PI, PI, 5);
    TF1 *fTimesCos2 = new TF1("fTimesCos2", AcceptanceFuncTimesCos, -PI, PI, 5);

    cout << "CORRECTIONS:\n";
    cout << "vn: {cm, sm, cm2, sm2, a-, a+, lambda-, lambda+}\n";

    double *corrections[nCoef];
    for (i=0; i<nCoef; i++) {
        corrections[i] = GetCorrectionParam(phiMin, phiMax, percentage, i+1, fA, fTimesSin, fTimesCos, fTimesSin2, fTimesCos2);
    }

    // Save correction parameters
    TH2D *hCorrectionParameters = new TH2D("hCorrectionParameters", "hCorrectionParameters", 5, 0.5, 5.5, 8, 0.5, 8.5);
    for (i=0; i<nCoef; i++) {
        for (j=0; j<nCorrParam; j++) {
            hCorrectionParameters->Fill(i, j, corrections[i][j]);
        }
    }

    // Save input numbers
    TH1D *hInputNumbers = new TH1D("hInputNumbers","hInputNumbers",14, 0.5, 14.5);
    hInputNumbers->Fill(1, double(nEvents));
    hInputNumbers->Fill(2, vn[0]);
    hInputNumbers->Fill(3, vn[1]);
    hInputNumbers->Fill(4, vn[2]);
    hInputNumbers->Fill(5, vn[3]);
    hInputNumbers->Fill(6, vn[4]);
    hInputNumbers->Fill(7, Tdec);
    hInputNumbers->Fill(8, vr);
    hInputNumbers->Fill(9, Teff);
    hInputNumbers->Fill(10, 1./Teff);
    hInputNumbers->Fill(11, phiMin);
    hInputNumbers->Fill(12, phiMax);
    hInputNumbers->Fill(13, percentage);
    hInputNumbers->Fill(14, 1.0); // Counting number of files added with hadd.
    fOut->cd();
    hInputNumbers->Write("hInputNumbers");

    int nOutput = nEvents/20;
    if (nOutput<1) nOutput = 1;
    for (i=0; i<nEvents; i++) {
        if (i % nOutput == 0)
            cout << 100*i/nEvents << " % finished" << endl;

        lists->ClearLists();

        if (bRandomPsi) {
            for (j=0; j<5; j++) {
                Psi[j] = rand->Uniform(-PI, PI);
            }
        } else {
            double psiTemp = rand->Uniform(-PI, PI);
            for (j=0; j<5;j++) {
                Psi[j] = psiTemp;
            }
        }

        fPhiDist->SetParameters(vn[0], vn[1], vn[2], vn[3], vn[4],
            Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

        GetEvent(histos, lists, inputs, rand, fPtDist, fPhiDist, fVnDist, vn, Psi, percentage, phiMin, phiMax, bUsePtDependence, bUseGranularity, -PI);
        AnalyzeEvent(lists, histos, Psi, bUseWeight, 0, corrections);
        AnalyzeEvent(lists, histos, Psi, bUseWeight, 1, corrections);

    }

    fOut->Write();
    fOut->Close();
    timer.Print();
    return 0;
}

//======END OF MAIN PROGRAM======
void GetEvent(JHistos *histos, JEventLists *lists, JInputs *inputs, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, TF1 *fVnDist, double *vn, double *Psi, double percentage, double phiMin, double phiMax, bool bUsePtDependence, bool bUseGranularity, double startAngle) {
    double pT, phi, eta, Energy;
    double px, py, pz;
    double centrality, randNum;
    double vnTemp[nCoef];

    double nMult, nMultNonuni;
    double alpha = 2.0, beta = 1.0;

    JToyMCTrack track;
    TLorentzVector lVec;
    TClonesArray *tempList;

    int i, j, k;

    centrality = rand->Uniform(0.0, 70.0);

    for (i=0; i<DET_N; i++) {

        nMult = 0;
        nMultNonuni = 0;

        for (j=0; j<CENTBINS_N; j++) {
            if (centrality>centBins[j] && centrality<centBins[j+1]) {
                centrality = centBins[j+1] - (centBins[j+1]-centBins[j])/2.0;
                histos->hCentrality->Fill(centrality);
                nMult = inputs->GetMultiplicity(i, centrality);
            }
        }

        for (j=0; j<nMult; j++) {
            eta = inputs->GetEta(i);
            pT = fPt->GetRandom();
            if (bUsePtDependence) {
                for (k=0; k<nCoef; k++) {
                    fVnDist->SetParameters(alpha, beta, vn[k]);
                    vnTemp[k] = fVnDist->Eval(pT);
                }
                fPhi->SetParameters(vnTemp[0], vnTemp[1], vnTemp[2], vnTemp[3], vnTemp[4],
                    Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);
            }

            phi = fPhi->GetRandom();
            if (bUseGranularity) {
                phi = CheckPhi(phi, startAngle);
                eta = CheckEta(eta);
            }

            histos->hPt->Fill(pT);
            histos->hPhi->Fill(phi);
            histos->hEta->Fill(eta);

            px = pT*TMath::Cos(phi);
            py = pT*TMath::Sin(phi);
            pz = pT*TMath::SinH(eta);
            Energy = TMath::Sqrt(pT*pT + pz*pz);
            lVec.SetPxPyPzE(px, py, pz, Energy);
            track.SetTrack(lVec);

            tempList = lists->GetList(i, 0);
            new((*tempList)[j]) JToyMCTrack(track);

            tempList = lists->GetList(i, 1);
            if (phi < phiMin || phi > phiMax) {

                histos->hPhiNonuni->Fill(phi);

                new((*tempList)[nMultNonuni]) JToyMCTrack(track);
                nMultNonuni++;

            } else {
                randNum = rand->Rndm();
                if ( randNum > percentage ) {

                    histos->hPhiNonuni->Fill(phi);

                    new((*tempList)[nMultNonuni]) JToyMCTrack(track);
                    nMultNonuni++;

                }
            }

        }
        histos->hMultiplicity[i]->Fill(nMult);
        histos->hMultiplicityNonuni[i]->Fill(nMultNonuni);
    }
}

void AnalyzeEvent(JEventLists *lists, JHistos *histos, double *Psi, bool bUseWeight, bool bDoCorrections, double **corrections) {

    int nMultTPC, nMultT0PA, nMultT0PC, nMultV0P;

    if (bDoCorrections) {
        nMultTPC = lists->TPClistNonuni->GetEntriesFast();
        nMultT0PA = lists->T0PAlistNonuni->GetEntriesFast();
        nMultT0PC = lists->T0PClistNonuni->GetEntriesFast();
        nMultV0P = lists->V0PlistNonuni->GetEntriesFast();
    } else {
        nMultTPC = lists->TPClist->GetEntriesFast();
        nMultT0PA = lists->T0PAlist->GetEntriesFast();
        nMultT0PC = lists->T0PClist->GetEntriesFast();
        nMultV0P = lists->V0Plist->GetEntriesFast();
    }

    int i, j, k, n, jMax;
    double w = 1.0;
    double phi, pt;
    double cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus;
    double EventPlaneA, EventPlaneB, Rtrue, Rsub, vobs;

    TComplex QvecTPC, QvecT0PA, QvecT0PC, QvecV0P;
    TComplex unitVec = TComplex(0, 0);
    TComplex autocorr = TComplex(0, 0);

    double normTPC, normT0PA, normT0PC, normV0P;

    double QnQnA, QnAQnB;

    vector<vector<TComplex>> pTBinsQ;
    pTBinsQ.resize(PTBINS_N);

    JToyMCTrack *track;

    for (i=0; i<nCoef; i++) {

        n = i+1;

        cm = corrections[i][0];
        sm = corrections[i][1];
        aMinus = corrections[i][4];
        aPlus = corrections[i][5];
        lambdaMinus = corrections[i][6];
        lambdaPlus = corrections[i][7];

        QvecTPC = TComplex(0, 0);
        QvecT0PA = TComplex(0, 0);
        QvecT0PC = TComplex(0, 0);
        QvecV0P = TComplex(0, 0);

        vobs = 0.0;

        normTPC = 0.0; normT0PA = 0.0; normT0PC = 0.0; normV0P = 0.0;

        // Construct Q-vectors for the detectors
        for (j=0; j<nMultTPC; j++) {
            track = (JToyMCTrack*)lists->TPClist->At(j);
            CalculateQvector(track, unitVec, QvecTPC, normTPC, bUseWeight, bDoCorrections, n, w, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
        }

        if (nMultT0PA<nMultT0PC) {
            jMax = nMultT0PA;
        } else {
            jMax = nMultT0PC;
        }

        for (j=0; j<jMax; j++) {
            track = (JToyMCTrack*)lists->T0PAlist->At(j);
            CalculateQvector(track, unitVec, QvecT0PA, normT0PA, bUseWeight, bDoCorrections, n, w, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);

            track = (JToyMCTrack*)lists->T0PClist->At(j);
            CalculateQvector(track, unitVec, QvecT0PC, normT0PC, bUseWeight, bDoCorrections, n, w, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
        }

        for (j=0; j<nMultV0P; j++) {
            track = (JToyMCTrack*)lists->V0Plist->At(j);
            CalculateQvector(track, unitVec, QvecV0P, normV0P, bUseWeight, bDoCorrections, n, w, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
        }

        // Calculate vobs from TPC events
        for (j=0; j<nMultTPC; j++) {

            track = (JToyMCTrack*)lists->TPClist->At(j);
            phi = track->GetPhi();
            pt = track->GetPt();

            if (bUseWeight) w = pt;

            autocorr = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));

            QvecTPC -= autocorr;
            vobs += GetVnObs(QvecTPC, phi, n);
            QvecTPC += autocorr;
        }

        vobs /= nMultTPC;

        // Resolution parameter calculations
        Rtrue = TMath::Cos(n*(GetEventPlane(QvecTPC, n) - Psi[i]));

        EventPlaneA = GetEventPlane(QvecT0PA, n);
        EventPlaneB = GetEventPlane(QvecT0PC, n);
        Rsub = TMath::Cos(n*(EventPlaneA - EventPlaneB));

        // ALTERNATIVE EVENT PLANE METHOD
        normTPC = TMath::Sqrt(normTPC);
        normT0PA = TMath::Sqrt(normT0PA);
        normT0PC = TMath::Sqrt(normT0PC);
        normV0P = TMath::Sqrt(normV0P);

        QvecTPC /= normTPC; QvecT0PA /= normT0PA; QvecT0PC /= normT0PC; QvecV0P /= normV0P;

        QnQnA = QvecV0P*TComplex::Conjugate(QvecT0PA);
        QnQnA /= TComplex::Abs(QvecT0PA);

        QnAQnB = QvecT0PA*TComplex::Conjugate(QvecT0PC);
        QnAQnB /= TComplex::Abs(QvecT0PA);
        QnAQnB /= TComplex::Abs(QvecT0PC);

        if (bDoCorrections) {
            histos->hVnObsCorrected[i]->Fill(vobs);
            histos->hRtrueCorrected[i]->Fill(Rtrue);
            histos->hRsubCorrected[i]->Fill(Rsub);
            histos->hQnQnAcorrected[i]->Fill(QnQnA);
            histos->hQnAQnBcorrected[i]->Fill(QnAQnB);
        } else {
            histos->hVnObs[i]->Fill(vobs);
            histos->hRtrue[i]->Fill(Rtrue);
            histos->hRsub[i]->Fill(Rsub);
            histos->hQnQnA[i]->Fill(QnQnA);
            histos->hQnAQnB[i]->Fill(QnAQnB);
        }

        // Divide into pT-bins
        if (n==2 && !bDoCorrections) {
            double weight = 0.0;
            double norms[PTBINS_N];
            for (j=0; j<PTBINS_N; j++) norms[j] = 0;

            for (j=0; j<nMultTPC; j++) {

                track = (JToyMCTrack*)lists->TPClist->At(j);
                phi = track->GetPhi();
                pt = track->GetPt();

                if (bUseWeight) w = pt;

                unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));

                for (k=0; k<PTBINS_N; k++) {
                    if ((pTBins[k] <= pt) && (pTBins[k+1] > pt)) {
                        pTBinsQ[k].push_back(unitVec);
                        norms[k] += w*w;
                    }
                }
            }

            double l;
            for (j=0; j<PTBINS_N; j++) {
                QvecTPC = TComplex(0, 0);
                l = pTBinsQ[j].size();
                for (k=0; k<l; k++) {
                    QvecTPC += pTBinsQ[j][k];
                }
                weight = TMath::Sqrt(norms[j]);
                if (weight!=0) QvecTPC /= weight;
                histos->hV2ComplexPart->Fill(QvecTPC.Im()*QvecT0PA.Im());
                QnQnA = QvecTPC*TComplex::Conjugate(QvecT0PA);
                QnQnA /= TComplex::Abs(QvecT0PA);
                histos->hPtBin[j]->Fill(QnQnA);
                histos->hSqrtSumWeightsPtBins[j]->Fill(weight);
                pTBinsQ[j].clear();
            }

        }
    }

    if (bDoCorrections) {
        histos->hSqrtSumWeightsTPCNonuni->Fill(normTPC);
        histos->hSqrtSumWeightsT0PANonuni->Fill(normT0PA);
        histos->hSqrtSumWeightsT0PCNonuni->Fill(normT0PC);
        histos->hSqrtSumWeightsV0PNonuni->Fill(normV0P);
    } else {
        histos->hSqrtSumWeightsTPC->Fill(normTPC);
        histos->hSqrtSumWeightsT0PA->Fill(normT0PA);
        histos->hSqrtSumWeightsT0PC->Fill(normT0PC);
        histos->hSqrtSumWeightsV0P->Fill(normV0P);
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
    double phi = x[0];
    double phiMin = p[0];
    double phiMax = p[1];
    double missedParticles = p[2];
    double norm = p[3];
    return ((phi < phiMin) || (phi > phiMax)) ? norm/(2*PI) : norm*(1-missedParticles)/(2*PI);
}

double AcceptanceFuncTimesSin(double *x, double *p) {
    double phi = x[0];
    double phiMin = p[0];
    double phiMax = p[1];
    double missedParticles = p[2];
    double n = p[3];
    double norm = p[4];
    return ((phi < phiMin) || (phi > phiMax)) ? norm*TMath::Sin(n*phi)/(2*PI) : norm*(1-missedParticles)*TMath::Sin(n*phi)/(2*PI);
}

double AcceptanceFuncTimesCos(double *x, double *p) {
    double phi = x[0];
    double phiMin = p[0];
    double phiMax = p[1];
    double missedParticles = p[2];
    double n = p[3];
    double norm = p[4];
    return ((phi < phiMin) || (phi > phiMax)) ? norm*TMath::Cos(n*phi)/(2*PI) : norm*(1-missedParticles)*TMath::Cos(n*phi)/(2*PI);
}

double *GetCorrectionParam(double phiMin, double phiMax, double percentage, double n, TF1 *fA, TF1 *fTimesSin, TF1 *fTimesCos, TF1 *fTimesSin2, TF1 *fTimesCos2) {
    fA->SetParameters(phiMin, phiMax, percentage, 1);
    double area = fA->Integral(0.0, 2*PI);

    fTimesSin->SetParameters(phiMin, phiMax, percentage, n, 1.0/area);
    fTimesCos->SetParameters(phiMin, phiMax, percentage, n, 1.0/area);
    fTimesSin2->SetParameters(phiMin, phiMax, percentage, 2*n, 1.0/area);
    fTimesCos2->SetParameters(phiMin, phiMax, percentage, 2*n, 1.0/area);

    double thres = 0.000001;
    double *corrections = new double[nCorrParam];
    corrections[0] = CheckIfZero(fTimesCos->Integral(-PI, PI), thres); //cm
    corrections[1] = CheckIfZero(fTimesSin->Integral(-PI, PI), thres); //sm
    corrections[2] = CheckIfZero(fTimesCos2->Integral(-PI, PI), thres); //cm2
    corrections[3] = CheckIfZero(fTimesSin2->Integral(-PI, PI), thres); //sm2
    corrections[4] = CheckIfZero(1.0 - corrections[2], thres); //a-
    corrections[5] = CheckIfZero(1.0 + corrections[2], thres); //a+
    corrections[6] = CheckIfZero(corrections[3]/corrections[4], thres); //lambda-
    corrections[7] = CheckIfZero(corrections[3]/corrections[5], thres); //lambda+

    cout << "v" << n << ": {";
    for (int i=0; i<nCorrParam-1; i++) {
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

void DoCorrections(TComplex &Qvec, double cm, double sm, double lambdaMinus, double lambdaPlus, double aMinus, double aPlus) {
    Qvec = TComplex(ShiftCorrection(Qvec.Re(), cm), ShiftCorrection(Qvec.Im(), sm));
    Qvec = TComplex(TwistCorrection(Qvec.Re(), Qvec.Im(), lambdaMinus, lambdaPlus), TwistCorrection(Qvec.Im(), Qvec.Re(), lambdaPlus, lambdaMinus));
    Qvec = TComplex(RescalingCorrection(Qvec.Re(), aPlus), RescalingCorrection(Qvec.Im(), aMinus));
}

void CalculateQvector(JToyMCTrack *track, TComplex unitVec, TComplex &Qvec, double &norm, bool bUseWeight, bool bDoCorrections, int n, double w, double cm, double sm, double lambdaMinus, double lambdaPlus, double aMinus, double aPlus) {
    double phi = track->GetPhi();
    double pt = track->GetPt();

    if (bUseWeight) w = pt;
    norm += w*w;

    unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
    if (bDoCorrections) DoCorrections(unitVec, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
    Qvec += unitVec;
}

double GetEventPlane(TComplex Qvec, int n) {
    return TMath::ATan2(Qvec.Im(), Qvec.Re())/n;
}

double GetVnObs(TComplex Qvec, double phi, int n) {
    return TMath::Cos(n*(phi - GetEventPlane(Qvec, n)));
}

double CheckPhi(double phi, double startAngle) {
    double lower, upper;
    for (int i=0; i<SECTORS_N; i++) {
        lower = startAngle+2*PI*i/SECTORS_N;
        upper = startAngle+2*PI*(i+1)/SECTORS_N;
        if (lower < phi && phi < upper) {
            return lower + (upper-lower)/2.0;
        }
    }
    return 0;
}

double CheckEta(double eta) {
    for (int i=0; i<RINGS_N; i++) {
        if (ringEta[i] < eta && eta < ringEta[i+1])
            return ringEta[i] + (ringEta[i+1]-ringEta[i])/2.0;
    }
    return 0;
}
