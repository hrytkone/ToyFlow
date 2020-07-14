#ifndef JEVENTLISTS_H
#define JEVENTLISTS_H

#include "TClonesArray.h"

class JEventLists {

public:
    JEventLists();
    virtual ~JEventLists() {;}

    void ClearLists();
    TClonesArray *GetList(int det_i, TString sAorB="");

    TClonesArray *fullEvent;
    TClonesArray *TPClist;
    TClonesArray *T0Alist;
    TClonesArray *T0Clist;
    TClonesArray *V0Alist;
    TClonesArray *TPCAlist;
    TClonesArray *TPCClist;
    TClonesArray *V0Clist;

    TClonesArray *TPClistA;
    TClonesArray *TPClistB;
    TClonesArray *T0AlistA;
    TClonesArray *T0AlistB;
    TClonesArray *T0ClistA;
    TClonesArray *T0ClistB;
    TClonesArray *V0AlistA;
    TClonesArray *V0AlistB;
    TClonesArray *TPCAlistA;
    TClonesArray *TPCAlistB;
    TClonesArray *TPCClistA;
    TClonesArray *TPCClistB;
    TClonesArray *V0ClistA;
    TClonesArray *V0ClistB;
};

#endif
