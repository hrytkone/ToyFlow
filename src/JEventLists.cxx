#include "JEventLists.h"
#include "JToyMCTrack.h"

JEventLists::JEventLists() {

    int size = 3000;

    fullEvent = new TClonesArray("JToyMCTrack", 10*size);
    TPClist = new TClonesArray("JToyMCTrack", size);
    T0Alist = new TClonesArray("JToyMCTrack", size);
    T0Clist = new TClonesArray("JToyMCTrack", size);
    V0Alist = new TClonesArray("JToyMCTrack", size);
    TPCAlist = new TClonesArray("JToyMCTrack", size);
    TPCClist = new TClonesArray("JToyMCTrack", size);
    V0Clist = new TClonesArray("JToyMCTrack", size);

    TPClistA = new TClonesArray("JToyMCTrack",size);
    TPClistB = new TClonesArray("JToyMCTrack",size);
    T0AlistA = new TClonesArray("JToyMCTrack",size);
    T0AlistB = new TClonesArray("JToyMCTrack",size);
    T0ClistA = new TClonesArray("JToyMCTrack",size);
    T0ClistB = new TClonesArray("JToyMCTrack",size);
    V0AlistA = new TClonesArray("JToyMCTrack",size);
    V0AlistB = new TClonesArray("JToyMCTrack",size);
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
    T0Alist->Clear("C");
    T0Clist->Clear("C");
    V0Alist->Clear("C");
    TPCAlist->Clear("C");
    TPCClist->Clear("C");
    V0Clist->Clear("C");
    TPClistA->Clear("C");
    TPClistB->Clear("C");
    T0AlistA->Clear("C");
    T0AlistB->Clear("C");
    T0ClistA->Clear("C");
    T0ClistB->Clear("C");
    V0AlistA->Clear("C");
    V0AlistB->Clear("C");
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
            case 1 : return T0AlistA;
            case 2 : return T0ClistA;
            case 3 : return V0AlistA;
            case 4 : return TPCAlistA;
            case 5 : return TPCClistA;
            case 6 : return V0ClistA;
        }
    } else if (sAorB.EqualTo("B")) {
        switch (det_i) {
            case 0 : return TPClistB;
            case 1 : return T0AlistB;
            case 2 : return T0ClistB;
            case 3 : return V0AlistB;
            case 4 : return TPCAlistB;
            case 5 : return TPCClistB;
            case 6 : return V0ClistB;
        }
    } else {
        switch (det_i) {
            case 0 : return TPClist;
            case 1 : return T0Alist;
            case 2 : return T0Clist;
            case 3 : return V0Alist;
            case 4 : return TPCAlist;
            case 5 : return TPCClist;
            case 6 : return V0Clist;
        }
    }
    return TPClist;
}
