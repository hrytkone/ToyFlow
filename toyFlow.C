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

void GetEvent(JHistos *histos, JEventLists *lists, JInputs *inputs, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, TF1 *fVnDist, double *vn, double *Psi, double percentage, double phiMin, double phiMax, bool bUsePtDependence, bool bUseGranularity, double startAngle, double centrality);
void GetParticleLists(JEventLists *lists, bool bDoCorrections);
void AnalyzeEvent(JHistos *histos, JEventLists *lists, JInputs *inputs, double *Psi, bool bUseWeight, bool bDoCorrections, double **corrections, double centrality);

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
double BelongsToA(double phi);

int main(int argc, char **argv) {

    TString outFileName = argc > 1 ? argv[1]:"toyFlow.root";
    if(outFileName.EqualTo("help",TString::kIgnoreCase)) {
        cout << "Usage: " << argv[0] << " filename.root nEvents bUsePtDep bUseGran seedNum" << endl;
        return 0;
    };
    int nEvents = argc > 2 ? atol(argv[2]):1000;
    bool bUsePtDependence = argc > 3 ? atol(argv[3]):0;
    bool bUseGranularity = argc > 4 ? atol(argv[4]):0;
    int iSeed = argc > 5 ? atol(argv[5]):0;

    bool bUseWeight = false;
    bool bRandomPsi = true;
    bool bUseCentDependence = false;

    const double scale = 1.0;
    double vn[nCoef] = {scale*0.0, scale*0.15, scale*0.08, scale*0.03, scale*0.01};

    cout << "=========================================== Settings ===========================================" << endl;
    cout << "Output: " << outFileName.Data()
         << ", Events: " << nEvents
         << ", Pt-dep: " << bUsePtDependence
         << ", Granularity: " << bUseGranularity
         << ", Seed: " << iSeed
         << ", Weight: " << bUseWeight
         << ", Random Psi: " << bRandomPsi
         << ", Centrality-dep: " << bUseCentDependence
         << endl;
    cout << "Vn inputs: ";
    for(int i=0; i<nCoef; i++) cout << vn[i] << ", ";
    cout << endl;
    cout << "================================================================================================" << endl;

    int i, j; // indices for loops

    TStopwatch timer;
    timer.Start();

    TFile *fOut = TFile::Open(outFileName, "RECREATE");

    gRandom->SetSeed(iSeed);
    TRandom3 *rand = new TRandom3(iSeed);
    //rand->SetSeed(0);

    JHistos *histos = new JHistos();
    JEventLists *lists = new JEventLists();
    JInputs *inputs = new JInputs();
    inputs->Load();

    double centrality;
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

    bool bDoCorrections;

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

        centrality = rand->Uniform(0.0, 60.0);

        if (bUseCentDependence) {
            for (j=0; j<nCoef; j++) vn[j] = inputs->GetCentDependVn(j+1, centrality);
        }

        fPhiDist->SetParameters(vn[0], vn[1], vn[2], vn[3], vn[4],
            Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

        GetEvent(histos, lists, inputs, rand, fPtDist, fPhiDist, fVnDist, vn, Psi, percentage, phiMin, phiMax, bUsePtDependence, bUseGranularity, -PI, centrality);

        // Analysis of uniform distribution
        bDoCorrections = false;
        GetParticleLists(lists, bDoCorrections);
        AnalyzeEvent(histos, lists, inputs, Psi, bUseWeight, bDoCorrections, corrections, centrality);

        // Analysis of non-uniform distribution
        bDoCorrections = true;
        GetParticleLists(lists, bDoCorrections);
        AnalyzeEvent(histos, lists, inputs, Psi, bUseWeight, bDoCorrections, corrections, centrality);

    }

    fOut->Write();
    fOut->Close();
    timer.Print();
    return 0;
}

//======END OF MAIN PROGRAM======
void GetEvent(JHistos *histos, JEventLists *lists, JInputs *inputs, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, TF1 *fVnDist, double *vn, double *Psi, double percentage, double phiMin, double phiMax, bool bUsePtDependence, bool bUseGranularity, double startAngle, double centrality) {
    double pT, phi, eta, Energy;
    double px, py, pz;
    double randNum;
    double vnTemp[nCoef];

    double nMult, nMultNonuni;
    double alpha = 2.0, beta = 1.0;

    JToyMCTrack track;
    TLorentzVector lVec;

    int centBin = 0;

    nMult = 0;
    nMultNonuni = 0;

    int i, j;
    for (i=0; i<CENTBINS_N-1; i++) {
        if (centrality>centBins[i] && centrality<centBins[i+1]) {
            centrality = centBins[i+1] - (centBins[i+1]-centBins[i])/2.0;
            centBin = i;
            histos->hCentrality->Fill(centrality);
            nMult = inputs->GetMultiplicity(centBin);
        }
    }

    for (i=0; i<nMult; i++) {
        eta = inputs->GetEta(centBin);
        pT = fPt->GetRandom();

        if (bUsePtDependence) {
            for (j=0; j<nCoef; j++) {
                fVnDist->SetParameters(alpha, beta, vn[j]);
                vnTemp[j] = fVnDist->Eval(pT);
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

        new((*lists->fullEvent)[i]) JToyMCTrack(track);

        if (phi < phiMin || phi > phiMax) {

            histos->hPhiNonuni->Fill(phi);

            new((*lists->fullEventNonuni)[nMultNonuni]) JToyMCTrack(track);
            nMultNonuni++;

        } else {
            randNum = rand->Rndm();
            if ( randNum > percentage ) {

                histos->hPhiNonuni->Fill(phi);

                new((*lists->fullEventNonuni)[nMultNonuni]) JToyMCTrack(track);
                nMultNonuni++;

            }
        }

    }

    histos->hMultiplicity->Fill(nMult);
    histos->hMultiplicityNonuni->Fill(nMultNonuni);

}

void GetParticleLists(JEventLists *lists, bool bDoCorrections) {
    int nMult = bDoCorrections ? lists->fullEventNonuni->GetEntriesFast()
                                : lists->fullEvent->GetEntriesFast();

    JToyMCTrack *tempTrack, track, trackA, trackB;
    double eta, phi;
    int detMult[DET_N] = {0};
    int detMultA[DET_N] = {0};
    int detMultB[DET_N] = {0};
    //double detCenter = 0.0;

    int i, j;
    for (i=0; i<nMult; i++) {
        tempTrack = bDoCorrections ? (JToyMCTrack*)lists->fullEventNonuni->At(i)
                                    : (JToyMCTrack*)lists->fullEvent->At(i);
        eta = tempTrack->GetEta();

        for (j=0; j<DET_N; j++) {
            if (cov[j][0]<eta && eta<cov[j][1]) {
                track = *tempTrack;
                new((*lists->GetList(j,bDoCorrections))[detMult[j]]) JToyMCTrack(track);
                detMult[j]++;

                phi = tempTrack->GetPhi();
                //detCenter = (cov[j][0]+cov[j][1])/2.0;
                if(/*eta<detCenter*/ i%2==0 /*BelongsToA(phi)*/) { //Later use the function.
                    trackA = *tempTrack;
                    new((*lists->GetList(j,bDoCorrections,"A"))[detMultA[j]]) JToyMCTrack(trackA);
                    detMultA[j]++;
                } else {
                    trackB = *tempTrack;
                    new((*lists->GetList(j,bDoCorrections,"B"))[detMultB[j]]) JToyMCTrack(trackB);
                    detMultB[j]++;
                }
            }
        }
    }
}

void AnalyzeEvent(JHistos *histos, JEventLists *lists, JInputs *inputs, double *Psi, bool bUseWeight, bool bDoCorrections, double **corrections, double centrality) {

    int nMult[DET_N], nMultA[DET_N], nMultB[DET_N];

    for(int iDet=0; iDet<DET_N; iDet++) {
        nMult[iDet] = lists->GetList(iDet,bDoCorrections)->GetEntriesFast();
        nMultA[iDet] = lists->GetList(iDet,bDoCorrections,"A")->GetEntriesFast();
        nMultB[iDet] = lists->GetList(iDet,bDoCorrections,"B")->GetEntriesFast();
    }

    int i, j, k, n, jMax;
    double w = 1.0;
    double phi, pt;
    double cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus;
    double EventPlaneA, EventPlaneB, Rtrue, Rsub, vobs;

    TComplex Qvec[DET_N], QvecA[DET_N], QvecB[DET_N];
    TComplex unitVec = TComplex(0, 0);
    TComplex autocorr = TComplex(0, 0);

    double norm[DET_N], normA[DET_N], normB[DET_N];

    double QnQnA, QnAQnB;

    int centBin = inputs->GetCentBin(centrality);

    vector<vector<TComplex>> pTBinsQ;
    pTBinsQ.resize(PTBINS_N);

    JToyMCTrack *track;;

    for (i=0; i<nCoef; i++) {

        n = i+1;
        cm = corrections[i][0];
        sm = corrections[i][1];
        aMinus = corrections[i][4];
        aPlus = corrections[i][5];
        lambdaMinus = corrections[i][6];
        lambdaPlus = corrections[i][7];

        for(int iDet=0; iDet<DET_N; iDet++) {

            Qvec[iDet]= TComplex(0, 0);
            QvecA[iDet] = TComplex(0, 0);
            QvecB[iDet] = TComplex(0, 0);

            vobs = 0.0;

            norm[iDet]= 0.0; normA[iDet] = 0.0; normB[iDet] = 0.0;

            // Construct Q-vectors for the detectors
            for (j=0; j<nMult[iDet]; j++) {
                track = (JToyMCTrack*)lists->GetList(iDet,bDoCorrections)->At(j);
                CalculateQvector(track, unitVec, Qvec[iDet], norm[iDet], bUseWeight, bDoCorrections, n, w, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
            }

            if (nMultA[iDet]<nMultB[iDet]) {
                jMax = nMultA[iDet];
            } else {
                jMax = nMultB[iDet];
            }

            for (j=0; j<jMax; j++) {
            //for (j=0; j<nMult[1]; j++) {
                //track = (JToyMCTrack*)lists->GetList(1,bDoCorrections)->At(j);
                track = (JToyMCTrack*)lists->GetList(iDet,bDoCorrections,"A")->At(j);
                CalculateQvector(track, unitVec, QvecA[iDet], normA[iDet], bUseWeight, bDoCorrections, n, w, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
            //}
            //for (j=0; j<nMult[2]; j++) {
                //track = (JToyMCTrack*)lists->GetList(2,bDoCorrections)->At(j);
                track = (JToyMCTrack*)lists->GetList(iDet,bDoCorrections,"B")->At(j);
                CalculateQvector(track, unitVec, QvecB[iDet], normB[iDet], bUseWeight, bDoCorrections, n, w, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
            }

            // Calculate vobs from TPC events
            for (j=0; j<nMult[iDet]; j++) {

                track = (JToyMCTrack*)lists->GetList(iDet,bDoCorrections)->At(j);
                phi = track->GetPhi();
                pt = track->GetPt();

                if (bUseWeight) w = pt;

                autocorr = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));

                Qvec[iDet] -= autocorr;
                vobs += GetVnObs(Qvec[iDet], phi, n);
                Qvec[iDet] += autocorr;
            }

            vobs /= nMult[iDet];

            // Resolution parameter calculations
            Rtrue = TMath::Cos(n*(GetEventPlane(Qvec[iDet], n) - Psi[i]));

            EventPlaneA = GetEventPlane(QvecA[iDet], n);
            EventPlaneB = GetEventPlane(QvecB[iDet], n);
            Rsub = TMath::Cos(n*(EventPlaneA - EventPlaneB));

            // EP-method and SP-method
            norm[iDet] = TMath::Sqrt(norm[iDet]);
            normA[iDet] = TMath::Sqrt(normA[iDet]);
            normB[iDet] = TMath::Sqrt(normB[iDet]);

            Qvec[iDet] /= norm[iDet]; QvecA[iDet] /= normA[iDet]; QvecB[iDet] /= normB[iDet];

            QnQnA = Qvec[iDet]*TComplex::Conjugate(QvecA[iDet]);
            //QnQnA /= TComplex::Abs(QvecA[iDet]);

            QnAQnB = QvecA[iDet]*TComplex::Conjugate(QvecB[iDet]);
            //QnAQnB /= TComplex::Abs(QvecA[iDet]);
            //QnAQnB /= TComplex::Abs(QvecB[iDet]);

            if (bDoCorrections) {
                histos->hVnObsNonuni[i][iDet][centBin]->Fill(vobs);
                histos->hRtrueNonuni[i][iDet][centBin]->Fill(Rtrue);
                histos->hRsubNonuni[i][iDet][centBin]->Fill(Rsub);
                histos->hQnQnAEPnonuni[i][iDet][centBin]->Fill(QnQnA/TComplex::Abs(QvecA[iDet]));
                histos->hQnAQnBEPnonuni[i][iDet][centBin]->Fill(QnAQnB/(TComplex::Abs(QvecA[iDet])*TComplex::Abs(QvecB[iDet])));
                histos->hQnQnASPnonuni[i][iDet][centBin]->Fill(QnQnA);
                histos->hQnAQnBSPnonuni[i][iDet][centBin]->Fill(QnAQnB);
            } else {
                histos->hVnObs[i][iDet][centBin]->Fill(vobs);
                histos->hRtrue[i][iDet][centBin]->Fill(Rtrue);
                histos->hRsub[i][iDet][centBin]->Fill(Rsub);
                histos->hQnQnAEP[i][iDet][centBin]->Fill(QnQnA/TComplex::Abs(QvecA[iDet]));
                histos->hQnAQnBEP[i][iDet][centBin]->Fill(QnAQnB/(TComplex::Abs(QvecA[iDet])*TComplex::Abs(QvecB[iDet])));
                histos->hQnQnASP[i][iDet][centBin]->Fill(QnQnA);
                histos->hQnAQnBSP[i][iDet][centBin]->Fill(QnAQnB);
            }
        }

        // Divide into pT-bins
        if (n==2 && !bDoCorrections) {
            double weight = 0.0;
            double norms[PTBINS_N];
            for (j=0; j<PTBINS_N; j++) norms[j] = 0;

            for (j=0; j<nMult[0]; j++) { //Only do this for TPC.

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
                Qvec[0] = TComplex(0, 0);
                l = pTBinsQ[j].size();
                for (k=0; k<l; k++) {
                    Qvec[0] += pTBinsQ[j][k];
                }
                weight = TMath::Sqrt(norms[j]);
                if (weight!=0) Qvec[0] /= weight;
                histos->hV2ComplexPart->Fill(Qvec[0].Im()*Qvec[1].Im());
                QnQnA = Qvec[0]*TComplex::Conjugate(Qvec[1]);
                QnQnA /= TComplex::Abs(Qvec[1]);
                histos->hQnQnAPtBin[j]->Fill(QnQnA);
                histos->hSqrtSumWeightsPtBins[j]->Fill(weight);
                pTBinsQ[j].clear();
            }

        }
    }

    for(int iDet=0; iDet<DET_N; iDet++) {
        if (bDoCorrections) {
            histos->hSqrtSumWeightsNonuni[iDet][centBin]->Fill(norm[iDet]);
        } else {
            histos->hSqrtSumWeights[iDet][centBin]->Fill(norm[iDet]);
        }
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
    return eta;
}

double BelongsToA(double phi) {
    return phi>0;
}
