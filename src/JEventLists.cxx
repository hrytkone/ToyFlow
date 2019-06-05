#include "JEventLists.h"
#include "JToyMCTrack.h"

JEventLists::JEventLists() {

    int size = 3000;

    fullEvent = new TClonesArray("JToyMCTrack", 10*size);
    TPClist = new TClonesArray("JToyMCTrack", size);
    T0PAlist = new TClonesArray("JToyMCTrack", size);
    T0PClist = new TClonesArray("JToyMCTrack", size);
    V0Plist = new TClonesArray("JToyMCTrack", size);

    fullEventNonuni = new TClonesArray("JToyMCTrack", 10*size);
    TPClistNonuni = new TClonesArray("JToyMCTrack", size);
    T0PAlistNonuni = new TClonesArray("JToyMCTrack", size);
    T0PClistNonuni = new TClonesArray("JToyMCTrack", size);
    V0PlistNonuni = new TClonesArray("JToyMCTrack", size);

    TPClistA = new TClonesArray("JToyMCTrack",size);
    TPClistB = new TClonesArray("JToyMCTrack",size);
    T0PAlistA = new TClonesArray("JToyMCTrack",size);
    T0PAlistB = new TClonesArray("JToyMCTrack",size);
    T0PClistA = new TClonesArray("JToyMCTrack",size);
    T0PClistB = new TClonesArray("JToyMCTrack",size);
    V0PlistA = new TClonesArray("JToyMCTrack",size);
    V0PlistB = new TClonesArray("JToyMCTrack",size);

    TPClistANonuni = new TClonesArray("JToyMCTrack",size);
    TPClistBNonuni = new TClonesArray("JToyMCTrack",size);
    T0PAlistANonuni = new TClonesArray("JToyMCTrack",size);
    T0PAlistBNonuni = new TClonesArray("JToyMCTrack",size);
    T0PClistANonuni = new TClonesArray("JToyMCTrack",size);
    T0PClistBNonuni = new TClonesArray("JToyMCTrack",size);
    V0PlistANonuni = new TClonesArray("JToyMCTrack",size);
    V0PlistBNonuni = new TClonesArray("JToyMCTrack",size);
}

void JEventLists::ClearLists() {
    fullEvent->Clear("C");
    TPClist->Clear("C");
    T0PAlist->Clear("C");
    T0PClist->Clear("C");
    V0Plist->Clear("C");
    fullEventNonuni->Clear("C");
    TPClistNonuni->Clear("C");
    T0PAlistNonuni->Clear("C");
    T0PClistNonuni->Clear("C");
    V0PlistNonuni->Clear("C");
    TPClistA->Clear("C");
    TPClistB->Clear("C");
    T0PAlistA->Clear("C");
    T0PAlistB->Clear("C");
    T0PClistA->Clear("C");
    T0PClistB->Clear("C");
    V0PlistA->Clear("C");
    V0PlistB->Clear("C");
    TPClistANonuni->Clear("C");
    TPClistBNonuni->Clear("C");
    T0PAlistANonuni->Clear("C");
    T0PAlistBNonuni->Clear("C");
    T0PClistANonuni->Clear("C");
    T0PClistBNonuni->Clear("C");
    V0PlistANonuni->Clear("C");
    V0PlistBNonuni->Clear("C");
}

TClonesArray *JEventLists::GetList(int det_i, bool bNonuni, TString sAorB) {
    if (!bNonuni) {
        if(sAorB.EqualTo("A")) {
            switch (det_i) {
                case 0 : return TPClistA;
                case 1 : return T0PAlistA;
                case 2 : return T0PClistA;
                case 3 : return V0PlistA;
            }
        } else if (sAorB.EqualTo("B")) {
            switch (det_i) {
                case 0 : return TPClistB;
                case 1 : return T0PAlistB;
                case 2 : return T0PClistB;
                case 3 : return V0PlistB;
            }
        } else {
            switch (det_i) {
                case 0 : return TPClist;
                case 1 : return T0PAlist;
                case 2 : return T0PClist;
                case 3 : return V0Plist;
            }
        }
    } else {
        if(sAorB.EqualTo("A")) {
            switch (det_i) {
                case 0 : return TPClistANonuni;
                case 1 : return T0PAlistANonuni;
                case 2 : return T0PClistANonuni;
                case 3 : return V0PlistANonuni;
            }
        } else if (sAorB.EqualTo("B")) {
            switch (det_i) {
                case 0 : return TPClistBNonuni;
                case 1 : return T0PAlistBNonuni;
                case 2 : return T0PClistBNonuni;
                case 3 : return V0PlistBNonuni;
            }
        } else {
            switch (det_i) {
                case 0 : return TPClistNonuni;
                case 1 : return T0PAlistNonuni;
                case 2 : return T0PClistNonuni;
                case 3 : return V0PlistNonuni;
            }
        }
    }
    return TPClist;
}
