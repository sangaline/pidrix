#ifndef Pidrixter_h
#define Pidrixter_h

#include "TNamed.h"
class Pidrix;

class Pidrixter : public TNamed {
  public:
    Pidrixter(const Pidrix *P = 0, unsigned int Nstat = 10, unsigned int Nsyst = 10);
    void Populate(const Pidrix *P, unsigned int Nstat = 10, unsigned int Nsyst = 10);
    void Reset();

    ClassDef(Pidrixter,1) //PID Matrix Factorization
};

#endif
