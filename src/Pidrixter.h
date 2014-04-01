#ifndef Pidrixter_h
#define Pidrixter_h

#include "TNamed.h"
class Pidrix;
class TObjArray;

class Pidrixter : public TNamed {
  private:
    TObjArray *pidrixes;
    unsigned int Nstat, Nsyst, Ntotal;

  public:
    Pidrixter(const Pidrix *P = 0, unsigned int Nstatistical = 10, unsigned int Nsystematic = 10);
    ~Pidrixter();
    void Populate(const Pidrix *P, unsigned int Nstatistical = 10, unsigned int Nsystematic = 10);
    void Reset();
    unsigned int Apply(bool (*updater)(Pidrix*));

    unsigned int Members() { return Ntotal; }
    Pidrix* Member(unsigned int index);


    ClassDef(Pidrixter,1) //PID Matrix Factorization
};

#endif
