#include "JEventLists.h"
#include "JToyMCTrack.h"

JEventLists::JEventLists() {

    int size = 3000;

    fullEvent = new TClonesArray("JToyMCTrack", 10*size);
    TPClist = new TClonesArray("JToyMCTrack", size);
    T0PAlist = new TClonesArray("JToyMCTrack", size);
    T0PClist = new TClonesArray("JToyMCTrack", size);
    V0Plist = new TClonesArray("JToyMCTrack", size);
    TPCAlist = new TClonesArray("JToyMCTrack", size);
    TPCClist = new TClonesArray("JToyMCTrack", size);
    V0Clist = new TClonesArray("JToyMCTrack", size);

    TPClistA = new TClonesArray("JToyMCTrack",size);
    TPClistB = new TClonesArray("JToyMCTrack",size);
    T0PAlistA = new TClonesArray("JToyMCTrack",size);
    T0PAlistB = new TClonesArray("JToyMCTrack",size);
    T0PClistA = new TClonesArray("JToyMCTrack",size);
    T0PClistB = new TClonesArray("JToyMCTrack",size);
    V0PlistA = new TClonesArray("JToyMCTrack",size);
    V0PlistB = new TClonesArray("JToyMCTrack",size);
    TPCAlistA = new TClonesArray("JToyMCTrack", size);
    TPCAlistB = new TClonesArray("JToyMCTrack", size);
    TPCClistA = new TClonesArray("JToyMCTrack", size);
    TPCClistB = new TClonesArray("JToyMCTrack", size);
    V0ClistA = new TClonesArray("JToyMCTrack", size);
    V0ClistB = new TClonesArray("JToyMCTrack", size);
}

void JEventLists::ClearLists() {
    fullEvent->Clear("C");
    TPClist->Clear("C");
    T0PAlist->Clear("C");
    T0PClist->Clear("C");
    V0Plist->Clear("C");
    TPCAlist->Clear("C");
    TPCClist->Clear("C");
    V0Clist->Clear("C");
    TPClistA->Clear("C");
    TPClistB->Clear("C");
    T0PAlistA->Clear("C");
    T0PAlistB->Clear("C");
    T0PClistA->Clear("C");
    T0PClistB->Clear("C");
    V0PlistA->Clear("C");
    V0PlistB->Clear("C");
    TPCAlistA->Clear("C");
    TPCAlistB->Clear("C");
    TPCClistA->Clear("C");
    TPCClistB->Clear("C");
    V0ClistA->Clear("C");
    V0ClistB->Clear("C");
}

TClonesArray *JEventLists::GetList(int det_i, TString sAorB) {
    if(sAorB.EqualTo("A")) {
        switch (det_i) {
            case 0 : return TPClistA;
            case 1 : return T0PAlistA;
            case 2 : return T0PClistA;
            case 3 : return V0PlistA;
            case 4 : return TPCAlistA;
            case 5 : return TPCClistA;
            case 6 : return V0ClistA;
        }
    } else if (sAorB.EqualTo("B")) {
        switch (det_i) {
            case 0 : return TPClistB;
            case 1 : return T0PAlistB;
            case 2 : return T0PClistB;
            case 3 : return V0PlistB;
            case 4 : return TPCAlistB;
            case 5 : return TPCClistB;
            case 6 : return V0ClistB;
        }
    } else {
        switch (det_i) {
            case 0 : return TPClist;
            case 1 : return T0PAlist;
            case 2 : return T0PClist;
            case 3 : return V0Plist;
            case 4 : return TPCAlist;
            case 5 : return TPCClist;
            case 6 : return V0Clist;
        }
    }
    return TPClist;
}
