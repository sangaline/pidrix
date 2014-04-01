#include "Pidrixter.h"
#include "Pidrix.h"

#include "TObjArray.h"

ClassImp(Pidrixter)

Pidrixter::Pidrixter(const Pidrix *P, unsigned int Nstatstical, unsigned int Nsystematic) {
    pidrixes = new TObjArray();
    pidrixes->SetOwner(kTRUE);

    if(P != 0) { 
        Populate(P, Nstatstical, Nsystematic);
    }
    else {
        Reset();
    }
}

Pidrixter::~Pidrixter() {
    delete pidrixes;
}

void Pidrixter::Reset() {
    Nstat = 0;
    Nsyst = 0;
    Ntotal = 0;
    pidrixes->Delete();
}

void Pidrixter::Populate(const Pidrix *P, unsigned int Nstatistical, unsigned int Nsystematic) {
    Reset();
    Nstat = Nstatistical;
    Nsyst = Nsystematic;
    Ntotal = Nstat*Nsyst;

    for(unsigned int i = 0; i < Nstat; i++) {
        //this Pidrix is a bootstrapped statistical variation on the main one
        Pidrix *P2 = new Pidrix();
        P2->SetTarget(P, true);
        P2->SetU(P);
        P2->SetV(P);
        P2->SetRank(P->Rank());
        pidrixes->AddLast(P2);
        for(unsigned int j = 1; j < Nsyst; j++) {
            //these are statistically identical but the initial conditions will vary
            Pidrix *P3 = new Pidrix();
            P3->SetTarget(P2);
            P3->SetU(P2);
            P3->SetV(P2);
            P3->SetRank(P->Rank());
            pidrixes->AddLast(P3);
        }
    }
}

unsigned int Pidrixter::Apply(bool (*updater)(Pidrix*)) {
    unsigned int positives = 0;
    for(unsigned int i = 0; i < Ntotal; i++) {
        Pidrix *P = (Pidrix*) pidrixes->At(i);
        if(updater(P)) positives++;
    }

    return positives;
}

Pidrix* Pidrixter::Member(unsigned int index) {
    return (Pidrix*) pidrixes->At(index);
}
