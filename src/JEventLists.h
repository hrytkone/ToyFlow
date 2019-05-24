#ifndef JEVENTLISTS_H
#define JEVENTLISTS_H

#include "TClonesArray.h"

class JEventLists {

public:
    JEventLists();
    virtual ~JEventLists() {;}

    void ClearLists();
    TClonesArray *GetList(int det_i, bool bNonuni);

    TClonesArray *TPClist;
    TClonesArray *T0PAlist;
    TClonesArray *T0PClist;
    TClonesArray *V0Plist;

    TClonesArray *TPClistNonuni;
    TClonesArray *T0PAlistNonuni;
    TClonesArray *T0PClistNonuni;
    TClonesArray *V0PlistNonuni;

};

#endif
