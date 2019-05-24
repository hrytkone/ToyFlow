#include "JEventLists.h"
#include "JToyMCTrack.h"

JEventLists::JEventLists() {

    int size = 3000;

    TPClist = new TClonesArray("JToyMCTrack", size);
    T0PAlist = new TClonesArray("JToyMCTrack", size);
    T0PClist = new TClonesArray("JToyMCTrack", size);
    V0Plist = new TClonesArray("JToyMCTrack", size);

    TPClistNonuni = new TClonesArray("JToyMCTrack", size);
    T0PAlistNonuni = new TClonesArray("JToyMCTrack", size);
    T0PClistNonuni = new TClonesArray("JToyMCTrack", size);
    V0PlistNonuni = new TClonesArray("JToyMCTrack", size);
}

void JEventLists::ClearLists() {
    TPClist->Clear("C");
    T0PAlist->Clear("C");
    T0PClist->Clear("C");
    V0Plist->Clear("C");
    TPClistNonuni->Clear("C");
    T0PAlistNonuni->Clear("C");
    T0PClistNonuni->Clear("C");
    V0PlistNonuni->Clear("C");
}

TClonesArray *JEventLists::GetList(int det_i, bool bNonuni) {
    if (!bNonuni) {
        switch (det_i) {
            case 0 : return TPClist;
            case 1 : return T0PAlist;
            case 2 : return T0PClist;
            case 3 : return V0Plist;
        }
    } else {
        switch (det_i) {
            case 0 : return TPClistNonuni;
            case 1 : return T0PAlistNonuni;
            case 2 : return T0PClistNonuni;
            case 3 : return V0PlistNonuni;
        }
    }
    return TPClist;
}
