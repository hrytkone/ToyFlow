#include "JEventLists.h"
#include "JToyMCTrack.h"

JEventLists::JEventLists() {

    int size = 3000;

    TPClist = new TClonesArray("JToyMCTrack", size);
    TPCAlist = new TClonesArray("JToyMCTrack", size);
    TPCClist = new TClonesArray("JToyMCTrack", size);
    V0Plist = new TClonesArray("JToyMCTrack", size);

    TPClistNonuni = new TClonesArray("JToyMCTrack", size);
    TPCAlistNonuni = new TClonesArray("JToyMCTrack", size);
    TPCClistNonuni = new TClonesArray("JToyMCTrack", size);
    V0PlistNonuni = new TClonesArray("JToyMCTrack", size);
}

void JEventLists::ClearLists() {
    TPClist->Clear("C");
    TPCAlist->Clear("C");
    TPCClist->Clear("C");
    V0Plist->Clear("C");
    TPClistNonuni->Clear("C");
    TPCAlistNonuni->Clear("C");
    TPCClistNonuni->Clear("C");
    V0PlistNonuni->Clear("C");
}

TClonesArray *JEventLists::GetList(int det_i, bool bNonuni) {
    if (!bNonuni) {
        switch (det_i) {
            case 0 : return TPClist;
            case 1 : return TPCAlist;
            case 2 : return TPCClist;
            case 3 : return V0Plist;
        }
    } else {
        switch (det_i) {
            case 0 : return TPClistNonuni;
            case 1 : return TPCAlistNonuni;
            case 2 : return TPCClistNonuni;
            case 3 : return V0PlistNonuni;
        }
    }
    return TPClist;
}
