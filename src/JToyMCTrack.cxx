#include "JToyMCTrack.h"
#include "TMath.h"

ClassImp(JToyMCTrack)

JToyMCTrack::JToyMCTrack() {
    lVec.SetXYZM(0.0, 0.0, 0.0, 0.0);
}

JToyMCTrack::JToyMCTrack(TLorentzVector lVec_in) {
    lVec.SetXYZM(lVec_in.Px(), lVec_in.Py(), lVec_in.Pz(), lVec_in.Mag());
}

JToyMCTrack::JToyMCTrack(TLorentzVector lVec_in, int ch_in, int pid_in, int isHadron_in) {
    lVec.SetXYZM(lVec_in.Px(), lVec_in.Py(), lVec_in.Pz(), lVec_in.Mag());
	lCharge = ch_in;
	lPID = pid_in;
	lIsHadron = isHadron_in;
}

void JToyMCTrack::SetMass(double m) {
    double px = lVec.Px();
    double py = lVec.Py();
    double pz = lVec.Pz();
    double E_in = TMath::Sqrt(m*m + px*px + py*py + pz*pz);
    lVec.SetE(E_in);
}

void JToyMCTrack::SetTrack(TLorentzVector lVec_in) {
    lVec = lVec_in;
}
