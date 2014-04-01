#include "Pidrix.h"

#include "plotting.h"

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TH2.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TRandom3.h"

ClassImp(Pidrix)

Pidrix::Pidrix() {
    T = new TMatrixD();
    U = new TMatrixD();
    V = new TMatrixD();
    E = new TMatrixD();
    SVD = new TDecompSVD();
    maxm = 0;
    maxn = 0;
    rank = 0;
    xlow = 0;
    xhigh = 0;
    ylow = 0;
    yhigh = 0;
    integral = 0;
}

Pidrix::~Pidrix() {
    delete T;
    delete U;
    delete V;
    delete E;
    delete SVD;
}

void Pidrix::SetTarget(const Pidrix *P, bool randomize) {
    TH2 *target = Plotting::Target(P,0);
    SetTarget(target, randomize);
    delete target;
}

void Pidrix::SetTarget(TH2 *target, bool randomize) {
    integral = target->Integral();
    if(randomize) {
        integral = gRandom->Poisson(integral);
        TH2* new_target = (TH2*) target->Clone();
        new_target->Reset();
        new_target->FillRandom(target, int(integral));
        target = new_target;
    }

    maxm = target->GetNbinsY();
    maxn = target->GetNbinsX();

    T->ResizeTo(maxm, maxn);
    E->ResizeTo(maxm, maxn);

    for(int i = 0; i < maxm; i++) {
        for(int j = 0; j < maxn; j++) {
            (*T)[i][j] = target->GetBinContent(j+1, i+1);
            (*E)[i][j] = target->GetBinError(j+1, i+1);
        }
    }
    xlow = target->GetXaxis()->GetBinLowEdge(1);
    xhigh = target->GetXaxis()->GetBinUpEdge(target->GetNbinsY());
    ylow = target->GetYaxis()->GetBinLowEdge(1);
    yhigh = target->GetYaxis()->GetBinUpEdge(target->GetNbinsY());

    //Prepare the SVD
    if(maxm < maxn) T->T();
    SVD->SetMatrix(*T);
    SVD->Decompose();
    if(maxm < maxn) T->T();

    //Reset the rank
    rank = 0;
    if(randomize) { delete target; }
}

unsigned int Pidrix::SetRank(unsigned int new_rank) { 
    U->ResizeTo(maxm, new_rank);
    V->ResizeTo(new_rank, maxn);
    
    return rank = new_rank; 
}

double Pidrix::SetU(unsigned i, unsigned j, double value) {
    return (*U)[i][j] = value;
}

double Pidrix::SetV(unsigned i, unsigned j, double value) {
    return (*V)[i][j] = value;
}
double Pidrix::GetU(unsigned i, unsigned j) {
    return (*U)[i][j];
}
double Pidrix::GetV(unsigned i, unsigned j) {
    return (*V)[i][j];
}
void Pidrix::SetU(const TMatrixD& newU) { 
    (*U) = newU; 
}
void Pidrix::SetV(const TMatrixD& newV) { 
    (*V) = newV;
}

void Pidrix::SetU(const Pidrix* P) {
    (*U) = P->GetU();
}

void Pidrix::SetV(const Pidrix* P) {
    (*V) = P->GetV();
}

const TVectorD& Pidrix::SVDSigma() const {
    return SVD->GetSig();
}

double Pidrix::RandomUniform(double min, double max) {
    return gRandom->Uniform(min, max);
}

