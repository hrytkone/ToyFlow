#ifndef JEVENTLISTS_H
#define JEVENTLISTS_H

#include "TClonesArray.h"

class JEventLists {

public:
    JEventLists();
    virtual ~JEventLists() {;}

    void ClearLists();
    TClonesArray *GetList(int det_i, bool bNonuni, TString sAorB="");

    TClonesArray *fullEvent;
    TClonesArray *TPClist;
    TClonesArray *T0PAlist;
    TClonesArray *T0PClist;
    TClonesArray *V0Plist;

    TClonesArray *fullEventNonuni;
    TClonesArray *TPClistNonuni;
    TClonesArray *T0PAlistNonuni;
    TClonesArray *T0PClistNonuni;
    TClonesArray *V0PlistNonuni;

    TClonesArray *TPClistA;
    TClonesArray *TPClistB;
    TClonesArray *T0PAlistA;
    TClonesArray *T0PAlistB;
    TClonesArray *T0PClistA;
    TClonesArray *T0PClistB;
    TClonesArray *V0PlistA;
    TClonesArray *V0PlistB;

    TClonesArray *TPClistANonuni;
    TClonesArray *TPClistBNonuni;
    TClonesArray *T0PAlistANonuni;
    TClonesArray *T0PAlistBNonuni;
    TClonesArray *T0PClistANonuni;
    TClonesArray *T0PClistBNonuni;
    TClonesArray *V0PlistANonuni;
    TClonesArray *V0PlistBNonuni;
};

#endif
