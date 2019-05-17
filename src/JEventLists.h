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
    TClonesArray *TPCAlist;
    TClonesArray *TPCClist;
    TClonesArray *V0Plist;

    TClonesArray *TPClistNonuni;
    TClonesArray *TPCAlistNonuni;
    TClonesArray *TPCClistNonuni;
    TClonesArray *V0PlistNonuni;

};

#endif
