#include "Pidrixter.h"
#include "Pidrix.h"

ClassImp(Pidrixter)

Pidrixter::Pidrixter(const Pidrix *P, unsigned int Nstat, unsigned int Nsyst) {
    if(P != 0) { 
        Populate(P, Nstat, Nsyst);
    }
}
