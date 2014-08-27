#include "Pidrixter.h"
#include "Pidrix.h"

#include "TObjArray.h"
#include "TMath.h"

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
        //this Pidrix is a bootstrapped statistical variation on the main one iff Nstat > 1
        Pidrix *P2 = new Pidrix();
        P2->SetTarget(P, Nstat > 1);
        P2->SetRank(P->Rank());
        P2->SetU(P);
        P2->SetV(P);
        pidrixes->AddLast(P2);
        for(unsigned int j = 1; j < Nsyst; j++) {
            //these are statistically identical but the initial conditions will vary
            Pidrix *P3 = new Pidrix();
            P3->SetTarget(P2);
            P3->SetRank(P->Rank());
            P3->SetU(P2);
            P3->SetV(P2);
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

double Pidrixter::Expectation(double (*eval)(Pidrix*)) {
    double sum = 0;
    for(unsigned int i = 0; i < Ntotal; i++) {
        sum += eval(Member(i));
    }
    return sum/double(Members());
}

double Pidrixter::TotalError(double (*eval)(Pidrix*)) {
    double sum = 0, sum2 = 0;
    for(unsigned int i = 0; i < Ntotal; i++) {
        const double value = eval(Member(i));
        sum += value;
        sum2 += value*value;
    }
    return TMath::Sqrt((sum2 - sum*sum/double(Ntotal))/double(Ntotal-1));
}

double Pidrixter::StatisticalError(double (*eval)(Pidrix*)) {
    return TMath::Sqrt( pow(TotalError(eval), 2) - pow(SystematicError(eval), 2) );
}

double Pidrixter::SystematicError(double (*eval)(Pidrix*)) {
    unsigned int index = 0;
    double average_error = 0;
    while(index < Ntotal) {
        double sum = 0, sum2 = 0;
        for(unsigned int j = 0; j < Nsyst; j++) {
            const double value = eval(Member(index));
            sum += value;
            sum2 += value*value;
            index++;
        }
        average_error += TMath::Sqrt((sum2 - sum*sum/double(Nsyst))/double(Nsyst-1));
    }
    return average_error/double(Nstat);
}
