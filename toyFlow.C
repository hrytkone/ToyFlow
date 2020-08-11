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
#include "TTree.h"
#include "TNtuple.h"

// OTHER
#include "TStopwatch.h"

using namespace std;

double PtDist(double *x, double *p);
double PhiDist(double *x, double *p);
double VnDist(double *x, double *p);

void GetEvent(JHistos *histos, JEventLists *lists, JInputs *inputs, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, TF1 *fVnDist, double *vn, double *Psi, double percentage, double phiMin, double phiMax, bool bNonuniformPhi, bool bUsePtDependence, double centrality, TNtuple *ntuple, int iEvt, double multiScale, double extraConvPart, double decays);
int AddParticle(JHistos *histos, JEventLists *lists, TRandom3 *rand, bool bNonuniformPhi, double phi, double phiMin, double phiMax, int listID, TLorentzVector lVec, int lCharge, int lPID, int lIsHadron, double percentage);
void GetParticleLists(JEventLists *lists, TRandom3* rand, bool bUseGranularity, bool bptCuts, double effTPC);
void AnalyzeEvent(JHistos *histos, JEventLists *lists, JInputs *inputs, double *Psi, bool bUseWeight, bool bNonuniformPhi, double **corrections, double centrality);
void AnalyzeUsing3sub(JHistos *histos, JEventLists *lists, JInputs *inputs, double centrality, bool bUseGranularity);

double AcceptanceFunc(double *x, double *p);
double AcceptanceFuncTimesSin(double *x, double *p);
double AcceptanceFuncTimesCos(double *x, double *p);
double *GetCorrectionParam(double phiMin, double phiMax, double percentage, double n, TF1 *fA, TF1 *fTimesSin, TF1 *fTimesCos, TF1 *fTimesSin2, TF1 *fTimesCos2);

double ShiftCorrection(double x, double correction);
double TwistCorrection(double x, double y, double lambda1, double lambda2);
double RescalingCorrection(double x, double a);
void DoCorrections(TComplex &Qvec, double cm, double sm, double lambdaMinus, double lambdaPlus, double aMinus, double aPlus);

void CalculateQvector(JToyMCTrack *track, TComplex unitVec, TComplex &Qvec, double &norm, bool bUseWeight, bool bNonuniformPhi, int n, double w, double cm, double sm, double lambdaMinus, double lambdaPlus, double aMinus, double aPlus);
double GetEventPlane(TComplex Qvec, int n);
double GetVnObs(TComplex Qvec, double phi, int n);

double CheckIfZero(double x, double thres);
double CheckPhi(double phi, double startAngle);
double CheckEta(double eta);
double BelongsToA(double phi);

// Arguments:
// filename.root filename of the output
// nEvents       number of events
// bUsePtDep     toggle flow pt dependency (currently not under development)
// scale         scale v_n signal strength
// multiScale    scale the multiplicity
// extraConvPart add a % of conversion particles (no flow)
// decays        decay a % of particles into two new particles (currently under construction)
// seedNum       seed number for the random number generator
// bSaveAsTrees  toggle between saving particles as trees with no histograms, or only histograms.
// bptCuts       do a pt cut at 150 MeV for TPC
// effTPC        set an efficiency for TPC so that particles at random are discarded.
int main(int argc, char **argv) {

    TString outFileName = argc > 1 ? argv[1]:"toyFlow.root";
    if(outFileName.EqualTo("help",TString::kIgnoreCase)) {
        cout << "Usage: " << argv[0] << " filename.root nEvents bUsePtDep bUseGran scale multiScale extraConvPart decays seedNum bSaveAsTrees bptCuts effTPC" << endl;
        return 0;
    };
    int nEvents = argc > 2 ? atol(argv[2]) : 1000;
    bool bUsePtDependence = argc > 3 ? atol(argv[3]) : 0;
    bool bUseGranularity = argc > 4 ? atol(argv[4]) : 0;
    double scale = argc > 5 ? atof(argv[5]) : 0.8;
    double multiScale = argc > 6 ? atof(argv[6]) : 1.0;
    double extraConvPart = argc > 7 ? atof(argv[7]) : 0.0;
    double decays = argc > 8 ? atof(argv[8]) : 0.0;
    int iSeed = argc > 9 ? atol(argv[9]) : 0;
    bool bSaveAsTrees = argc > 10 ? atol(argv[10]) : 0;
    bool bptCuts = argc > 11 ? atol(argv[11]) : 0;
    double effTPC = argc > 12 ? atof(argv[12]) : 1.0;

    bool bUseWeight = false;
    bool bRandomPsi = true;
    bool bUseCentDependence = true;
    bool bNonuniformPhi = false;

    double vn[nCoef] = {scale*0.0, scale*0.15, scale*0.08, scale*0.03, scale*0.01};

    cout << "=========================================== Settings ===========================================" << endl;
    cout << "Output: " << outFileName.Data()
         << ", Events: " << nEvents
         << ", Pt-dep: " << bUsePtDependence
         << ", Granularity: " << bUseGranularity
         << ", vn scale: " << scale
         << ", multiplicity scale: " << multiScale
         << ", extra conversion particles: " << extraConvPart
         << ", decay perc of the particles: " << decays
         << ", Seed: " << iSeed
         << ", Weight: " << bUseWeight
         << ", Random Psi: " << bRandomPsi
         << ", Centrality-dep: " << bUseCentDependence
         << ", Nonuniform phi: " << bNonuniformPhi
         << endl;
    cout << "Vn inputs: ";
    for(int i=0; i<nCoef; i++) cout << vn[i] << ", ";
    cout << endl;
    cout << "================================================================================================" << endl;

    int i, j; // indices for loops

    TStopwatch timer;
    timer.Start();

    JEventLists *lists = new JEventLists();
    JInputs *inputs = new JInputs();
    inputs->Load();

    TFile *fOut = TFile::Open(outFileName, "RECREATE");

    gRandom->SetSeed(iSeed);
    TRandom3 *rand = new TRandom3(iSeed);
    //rand->SetSeed(0);

    JHistos *histos = 0;
    TNtuple *ntuple = 0;
    if(bSaveAsTrees) {
        ntuple = new TNtuple("events", "data from toyFlow", "eventId:particleId:px:py:pz:x:y:z:isHadron");
    } else {
        // Do not produce histos if saving as trees.
        histos = new JHistos();
    }

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

        centrality = rand->Uniform(0.0, 50.0);

        if (bUseCentDependence) {
            for (j=0; j<nCoef; j++) vn[j] = scale * inputs->GetCentDependVn(j+1, centrality);
        }

        fPhiDist->SetParameters(vn[0], vn[1], vn[2], vn[3], vn[4],
            Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);

        GetEvent(histos, lists, inputs, rand, fPtDist, fPhiDist, fVnDist, vn, Psi, percentage, phiMin, phiMax, bNonuniformPhi, bUsePtDependence, centrality, ntuple, i, multiScale, extraConvPart, decays);
        if(bSaveAsTrees) {
            //
        } else { //When saving tracks as trees, we don't need to analyze event.
            GetParticleLists(lists, rand, bUseGranularity, bptCuts, effTPC);
            AnalyzeEvent(histos, lists, inputs, Psi, bUseWeight, bNonuniformPhi, corrections, centrality);
            AnalyzeUsing3sub(histos, lists, inputs, centrality, bUseGranularity);
        }
    }

    //Using the kOverwrite option a previous key with the same name is overwritten. The previous key is deleted before writing the new object.
    fOut->Write("",TObject::kOverwrite);
    fOut->Close();
    timer.Print();
    return 0;
}

//======END OF MAIN PROGRAM======
void GetEvent(JHistos *histos, JEventLists *lists, JInputs *inputs, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, TF1 *fVnDist, double *vn, double *Psi, double percentage, double phiMin, double phiMax, bool bNonuniformPhi, bool bUsePtDependence, double centrality, TNtuple *ntuple, int iEvt, double multiScale, double extraConvPart, double decays) {
    double pT, phi, eta, Energy;
    double pionMass = 0.135;
    double px, py, pz;
    double randNum;
    double vnTemp[nCoef];
    int lPID;
    int lCharge;
    int lIsHadron;

    int nMult = 0, nMultOrig =0, nTracks = 0;
    double alpha = 2.0, beta = 1.0;
    //http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf PID
    //Thermal model AN:
    //https://alice-notes.web.cern.ch/system/files/notes/analysis/308/2014-May-13-analysis_note-ThermalFitsAnalysisNote.pdf
    //check also source [5]: https://inspirehep.net/record/1222333
    // Ratios from 10-20 % centrality class.
    // pi+    pi-    K+   K-   p        anti-p
    // 455±31 453±31 68±5 68±6 21.0±1.7 21.1±1.8
    double probPiPlus  = 455.0/(455+455+68+68+21+21);
    double probPiMinus = 455.0/(455+455+68+68+21+21);
    double probKPlus   = 68.0/(455+455+68+68+21+21);
    double probKMinus  = 68.0/(455+455+68+68+21+21);
    double probP       = 21.0/(455+455+68+68+21+21);
    double probAntiP   = 21.0/(455+455+68+68+21+21);

    JToyMCTrack track;
    TLorentzVector lVec, lVecDecay1, lVecDecay2;
    TVector3 lBoost3Vec;

    int i, j, centBin = 0;
    for (i=0; i<CENTBINS_N-1; i++) {
        if (centrality>centBins[i] && centrality<centBins[i+1]) {
            centrality = centBins[i+1] - (centBins[i+1]-centBins[i])/2.0;
            centBin = i;
            if(histos!=0) histos->hCentrality->Fill(centrality);
            // Extra multiplicity scaling * extra conversion particles (no flow) * Multiplicity 
            nMultOrig = multiScale*inputs->GetMultiplicity(centBin);
            nMult = (1.0+extraConvPart)*nMultOrig;
        }
    }

    // Note that this i is incremented also in the decay if the particle happens to decay.
    for (i=0; i<nMult; i++) {

        eta = inputs->GetEta(centrality);
        pT = fPt->GetRandom(); // from the analytic function
        //pT = inputs->GetPt(centrality); // from pT centrality dependent distribution

        if (bUsePtDependence) {
            for (j=0; j<nCoef; j++) {
                fVnDist->SetParameters(alpha, beta, vn[j]);
                vnTemp[j] = fVnDist->Eval(pT);
            }
            fPhi->SetParameters(vnTemp[0], vnTemp[1], vnTemp[2], vnTemp[3], vnTemp[4],
                Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);
        }

        /* // 
        if ((eta>cov[D_TPC][0]) && (eta<cov[D_TPC][1])) {
             if (pT<0.2) continue;
             if (rand->Uniform()>0.8) continue;
        }
        */

        // The last extraConvPart % of particles are conversion particles with no flow.
        if(i > nMultOrig) {
            phi = rand->Uniform(-PI,PI);
        } else {
            phi = fPhi->GetRandom();
        }
        px = pT*TMath::Cos(phi);
        py = pT*TMath::Sin(phi);
        pz = pT*TMath::SinH(eta);
        Energy = TMath::Sqrt(pT*pT + pz*pz + pionMass*pionMass);
        lVec.SetPxPyPzE(px, py, pz, Energy);
        track.SetTrack(lVec);

        // Selecting particle species.
        randNum = rand->Rndm();
        if(randNum<probPiPlus) {
            lPID = 211; lCharge = 1; lIsHadron = 1; //pi+
        } else if(randNum<probPiPlus+probPiMinus) {
            lPID = -211; lCharge = -1; lIsHadron = 1; //pi-
        } else if(randNum<probPiPlus+probPiMinus+probKPlus) {
            lPID = 321; lCharge = 1; lIsHadron = 1; //K+
        } else if(randNum<probPiPlus+probPiMinus+probKPlus+probKMinus) {
            lPID = -321; lCharge = -1; lIsHadron = 1; //K-
        } else if(randNum<probPiPlus+probPiMinus+probKPlus+probKMinus+probP) {
            lPID = 2212; lCharge = 1; lIsHadron = 1; //p
        } else if(randNum<probPiPlus+probPiMinus+probKPlus+probKMinus+probP+probAntiP) {
            lPID = -2212; lCharge = -1; lIsHadron = 1; //anti-p
        } else {
            lPID = 0; lCharge = -9; lIsHadron = -1;
        }

        // When decays is over 0, that perc of the particles are decayed,
        // but the multiplicity is kept the same. This is for studying
        // the effects of strongly correlated decay particles in flow.
        randNum = rand->Rndm();
        if(randNum < decays) { // This needs to be checked randomly.
            // Calculate lorentz vectors for decay particles in CMS frame
            double decayEnergy = pionMass/2.0; // Pion mass halved -> assume pion mass for decaying particle.
            double phiDecay = rand->Uniform(0.0, 2*PI);
            double thetaDecay = rand->Uniform(0.0, PI);
            double pxDec = decayEnergy*TMath::Cos(phiDecay)*TMath::Sin(thetaDecay);
            double pyDec = decayEnergy*TMath::Sin(phiDecay)*TMath::Sin(thetaDecay);
            double pzDec = decayEnergy*TMath::Cos(thetaDecay);

            // decay products back to back in CMS frame
            lVecDecay1.SetPxPyPzE(pxDec, pyDec, pzDec, decayEnergy);
            lVecDecay2.SetPxPyPzE(-pxDec, -pyDec, -pzDec, decayEnergy);

            // Boost decay products to LAB frame
            lBoost3Vec = lVec.BoostVector();
            lVecDecay1.Boost(lBoost3Vec);
            lVecDecay2.Boost(lBoost3Vec);

            // Note that this decay is not taking into account any charge conservation or other such details.
            nTracks += AddParticle(histos, lists, rand, bNonuniformPhi, lVecDecay1.Phi(), phiMin, phiMax, nTracks, lVecDecay1, lCharge, lPID, lIsHadron, percentage);
            nTracks += AddParticle(histos, lists, rand, bNonuniformPhi, lVecDecay2.Phi(), phiMin, phiMax, nTracks, lVecDecay2, lCharge, lPID, lIsHadron, percentage);

            if(histos!=0) {
                histos->hPt->Fill(lVecDecay1.Pt());
                histos->hPt->Fill(lVecDecay2.Pt());
                histos->hEta->Fill(lVecDecay1.Eta());
                histos->hEta->Fill(lVecDecay2.Eta());
            }

            if(ntuple!=0) {
                ntuple->Fill(iEvt, lPID, lVecDecay1.Px(), lVecDecay1.Py(), lVecDecay1.Pz(), 0, 0, 0, lIsHadron);
                ntuple->Fill(iEvt, lPID, lVecDecay2.Px(), lVecDecay2.Py(), lVecDecay2.Pz(), 0, 0, 0, lIsHadron);
            }

            // One extra i++ is needed here as the decay results in two particles instead of the normal one particle.
            i++;
        } else {
            nTracks += AddParticle(histos, lists, rand, bNonuniformPhi, phi, phiMin, phiMax, nTracks, lVec, lCharge, lPID, lIsHadron, percentage);

            if(histos!=0) histos->hPt->Fill(pT);
            if(histos!=0) histos->hEta->Fill(eta);

            if(ntuple!=0) {
                ntuple->Fill(iEvt, lPID, px, py, pz, 0, 0, 0, lIsHadron);
            }
        }
    }

    if(histos!=0) histos->hMultiplicity[centBin]->Fill(nTracks);
}

// Returns 1 or 0 depending on if a particle was added to the list or not.
int AddParticle(JHistos *histos, JEventLists *lists, TRandom3 *rand, bool bNonuniformPhi, double phi, double phiMin, double phiMax, int listID, TLorentzVector lVec, int lCharge, int lPID, int lIsHadron, double percentage) {
    double randNum;
    int particleAdded = 0;
    if (phi < phiMin || phi > phiMax) {
        if(histos!=0) histos->hPhi->Fill(phi);
        // track, charge, pid, ishadron
        new((*lists->fullEvent)[listID]) JToyMCTrack(lVec, lCharge, lPID, lIsHadron);
        particleAdded = 1;
    } else {
        randNum = bNonuniformPhi ? rand->Rndm():1.1;
        if ( randNum > percentage ) {
            if(histos!=0) histos->hPhi->Fill(phi);
            // track, charge, pid, ishadron
            new((*lists->fullEvent)[listID]) JToyMCTrack(lVec, lCharge, lPID, lIsHadron);
            particleAdded = 1;
        }
    }
    return particleAdded;
}

void GetParticleLists(JEventLists *lists, TRandom3* rand, bool bUseGranularity, bool bptCuts, double effTPC) {
    int nMult = lists->fullEvent->GetEntriesFast();

    JToyMCTrack *tempTrack, track, trackA, trackB;
    double eta, phi, drand;
    bool discardParticleTPC = false;
    bool isTPC = false;
    int detMult[DET_N] = {0};
    int detMultA[DET_N] = {0};
    int detMultB[DET_N] = {0};
    //double detCenter = 0.0;

    int i, j;
    for (i=0; i<nMult; i++) {
        tempTrack = (JToyMCTrack*)lists->fullEvent->At(i);
        eta = tempTrack->GetEta();

        if(effTPC < 1.0) {
            drand = rand->Uniform();
            if(drand > effTPC) discardParticleTPC = true;
            else discardParticleTPC = false;
        }

        for (j=0; j<DET_N; j++) {

            isTPC = (j==D_TPC || j==D_TPC_A || j==D_TPC_C);
            if(isTPC && discardParticleTPC) continue;

            if (cov[j][0]<eta && eta<cov[j][1]) {
                phi = tempTrack->GetPhi();
                // Don't include tracks with less than 150 MeV in TPC.
                if(bptCuts && isTPC && tempTrack->GetPt() < 0.15 ) continue; 
                if (bUseGranularity && j==3) {
                    phi = CheckPhi(phi, -PI);
                    tempTrack->SetPhi(phi);
                    eta = CheckEta(eta); // need to add SetEta to JToyMCTrack.h
                }

                // Phi granularity for V0C
                if (bUseGranularity && j==6) {
                    phi = CheckPhi(phi, -PI);
                    tempTrack->SetPhi(phi);
                }

                track = *tempTrack;
                new((*lists->GetList(j))[detMult[j]]) JToyMCTrack(track);
                detMult[j]++;

                //detCenter = (cov[j][0]+cov[j][1])/2.0;
                if(/*eta<detCenter*/ i%2==0 /*BelongsToA(phi)*/) { //Later use the function.
                    trackA = *tempTrack;
                    new((*lists->GetList(j,"A"))[detMultA[j]]) JToyMCTrack(trackA);
                    detMultA[j]++;
                } else {
                    trackB = *tempTrack;
                    new((*lists->GetList(j,"B"))[detMultB[j]]) JToyMCTrack(trackB);
                    detMultB[j]++;
                }
            }
        }
    }
}

void AnalyzeEvent(JHistos *histos, JEventLists *lists, JInputs *inputs, double *Psi, bool bUseWeight, bool bNonuniformPhi, double **corrections, double centrality) {

    int nMult[DET_N], nMultA[DET_N], nMultB[DET_N];

    for(int iDet=0; iDet<DET_N; iDet++) {
        nMult[iDet] = lists->GetList(iDet)->GetEntriesFast();
        nMultA[iDet] = lists->GetList(iDet,"A")->GetEntriesFast();
        nMultB[iDet] = lists->GetList(iDet,"B")->GetEntriesFast();
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

    for(int iDet=0; iDet<DET_N; iDet++) {
        histos->hMultiPerDet[iDet][centBin]->Fill(nMult[iDet]);
    }

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
                track = (JToyMCTrack*)lists->GetList(iDet)->At(j);
                CalculateQvector(track, unitVec, Qvec[iDet], norm[iDet], bUseWeight, bNonuniformPhi, n, w, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
            }

            if (nMultA[iDet]<nMultB[iDet]) {
                jMax = nMultA[iDet];
            } else {
                jMax = nMultB[iDet];
            }

            for (j=0; j<jMax; j++) {
            //for (j=0; j<nMult[1]; j++) {
                //track = (JToyMCTrack*)lists->GetList(1)->At(j);
                track = (JToyMCTrack*)lists->GetList(iDet,"A")->At(j);
                CalculateQvector(track, unitVec, QvecA[iDet], normA[iDet], bUseWeight, bNonuniformPhi, n, w, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
            //}
            //for (j=0; j<nMult[2]; j++) {
                //track = (JToyMCTrack*)lists->GetList(2)->At(j);
                track = (JToyMCTrack*)lists->GetList(iDet,"B")->At(j);
                CalculateQvector(track, unitVec, QvecB[iDet], normB[iDet], bUseWeight, bNonuniformPhi, n, w, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
            }

            // Calculate vobs from TPC events! This is always done from TPC for each detector.
            for (j=0; j<nMult[D_TPC]; j++) {

                track = (JToyMCTrack*)lists->GetList(D_TPC)->At(j);
                phi = track->GetPhi();
                pt = track->GetPt();

                if (bUseWeight) w = pt;

                // No need to remove autocorrelation except for TPC as there the particles
                // which are used to calculate Qvec are also used to calculate phi.
                switch (iDet) {
                    case D_TPC : // Whole TPC needs to remove autocorrelation every time.
                        autocorr = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
                        break;
                    case D_TPC_A : // TPC_A needs to remove autocorrelation only if particle
                             // is inside TPC_A acceptance
                        if(track->GetEta()>cov[D_TPC_A][0] && track->GetEta()<cov[D_TPC_A][0]) {
                            autocorr = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
                        } else {
                            autocorr = TComplex(0,0);
                        }
                        break;
                    case D_TPC_C : // TPC_C needs to remove autocorrelation only if particle
                             // is inside TPC_C acceptance
                        if(track->GetEta()>cov[D_TPC_C][0] && track->GetEta()<cov[D_TPC_C][0]) {
                            autocorr = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
                        } else {
                            autocorr = TComplex(0,0);
                        }
                        break;
                    default : 
                        autocorr = TComplex(0,0);
                }

                Qvec[iDet] -= autocorr;
                vobs += GetVnObs(Qvec[iDet], phi, n);
                Qvec[iDet] += autocorr;
            }

            // TPC multi used to scale vobs
            vobs /= nMult[D_TPC];

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
            QnAQnB = QvecA[iDet]*TComplex::Conjugate(QvecB[iDet]);

            histos->hVnObs[i][iDet][centBin]->Fill(vobs);
            histos->hRtrue[i][iDet][centBin]->Fill(Rtrue);
            histos->hRsub[i][iDet][centBin]->Fill(Rsub);
            histos->hQnQnAEP[i][iDet][centBin]->Fill(QnQnA/TComplex::Abs(QvecA[iDet]));
            histos->hQnAQnBEP[i][iDet][centBin]->Fill(QnAQnB/(TComplex::Abs(QvecA[iDet])*TComplex::Abs(QvecB[iDet])));
            histos->hQnQnASP[i][iDet][centBin]->Fill(QnQnA);
            histos->hQnAQnBSP[i][iDet][centBin]->Fill(QnAQnB);
        }

        // Divide into pT-bins
        if (n==2) {
            double weight = 0.0;
            double norms[PTBINS_N];
            for (j=0; j<PTBINS_N; j++) norms[j] = 0;

            for (j=0; j<nMult[D_TPC]; j++) { //Only do this for TPC.

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
                Qvec[D_TPC] = TComplex(0, 0);
                l = pTBinsQ[j].size();
                for (k=0; k<l; k++) {
                    Qvec[D_TPC] += pTBinsQ[j][k];
                }
                weight = TMath::Sqrt(norms[j]);
                if (weight!=0) Qvec[D_TPC] /= weight;
                histos->hV2ComplexPart->Fill(Qvec[D_TPC].Im()*Qvec[1].Im());
                QnQnA = Qvec[D_TPC]*TComplex::Conjugate(Qvec[1]);
                QnQnA /= TComplex::Abs(Qvec[1]);
                histos->hQnQnAPtBin[j]->Fill(QnQnA);
                histos->hSqrtSumWeightsPtBins[j]->Fill(weight);
                pTBinsQ[j].clear();
            }

        }
    }

    for(int iDet=0; iDet<DET_N; iDet++)
        histos->hSqrtSumWeights[iDet][centBin]->Fill(norm[iDet]);
}

void AnalyzeUsing3sub(JHistos *histos, JEventLists *lists, JInputs *inputs, double centrality, bool bUseGranularity) {
    // 0: D_TPC   -0.8,  0.8
    // 1: D_T0_A   4.5,  5.0
    // 2: D_T0_C  -3.3, -2.9
    // 3: D_V0_A   2.8,  5.1
    // 4: D_TPC_A  0.1,  0.8
    // 5: D_TPC_C -0.8, -0.1
    // 6: D_V0_C  -3.7, -1.7
    const int detA = 6;
    const int detB = 4;
    const int detC = 5;

    TComplex QvecA, QvecB, QvecC;
    TComplex unitVec = TComplex(0, 0);

    double normA, normB, normC;

    JToyMCTrack *track;

    vector<double> phiList;

    int nmultA = lists->GetList(detA)->GetEntriesFast();
    int nmultB = lists->GetList(detB)->GetEntriesFast();
    int nmultC = lists->GetList(detC)->GetEntriesFast();

    for (Int_t j=0; j<nCoef; j++) {

        int n = j+1;

        QvecA = TComplex(0, 0);
        QvecB = TComplex(0, 0);
        QvecC = TComplex(0, 0);

        normA = 0.0; normB = 0.0; normC = 0.0;

        for (Int_t i=0; i<nmultA; i++) {
            track = (JToyMCTrack*)lists->GetList(detA)->At(i);
            CalculateQvector(track, unitVec, QvecA, normA, 0, 0, n, 1.0, 0, 0, 0, 0, 0, 0);
        }
        for (Int_t i=0; i<nmultB; i++) {
            track = (JToyMCTrack*)lists->GetList(detB)->At(i);
            CalculateQvector(track, unitVec, QvecB, normB, 0, 0, n, 1.0, 0, 0, 0, 0, 0, 0);
        }
        for (Int_t i=0; i<nmultC; i++) {
            track = (JToyMCTrack*)lists->GetList(detC)->At(i);
            CalculateQvector(track, unitVec, QvecC, normC, 0, 0, n, 1.0, 0, 0, 0, 0, 0, 0);
        }

        int ibin = inputs->GetCentBin(centrality);

        double epA = GetEventPlane(QvecA, n);
        double epB = GetEventPlane(QvecB, n);
        double epC = GetEventPlane(QvecC, n);
        histos->hRsubAB[j][ibin]->Fill(TMath::Cos((n)*(epA - epB)));
        histos->hRsubAC[j][ibin]->Fill(TMath::Cos((n)*(epA - epC)));
        histos->hRsubBC[j][ibin]->Fill(TMath::Cos((n)*(epB - epC)));
    }
}

double PtDist(double *x, double *p) {
    return x[0]*TMath::Exp(-p[0]*x[0]);
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

void CalculateQvector(JToyMCTrack *track, TComplex unitVec, TComplex &Qvec, double &norm, bool bUseWeight, bool bNonuniformPhi, int n, double w, double cm, double sm, double lambdaMinus, double lambdaPlus, double aMinus, double aPlus) {

    double phi = track->GetPhi();
    double pt = track->GetPt();

    if (bUseWeight) w = pt;
    norm += w*w;

    unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
    if (bNonuniformPhi) DoCorrections(unitVec, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
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
