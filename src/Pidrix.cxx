#include "Pidrix.h"

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TH2.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TRandom3.h"

ClassImp(Pidrix);

Pidrix::Pidrix() {
    T = new TMatrixD();
    U = new TMatrixD();
    V = new TMatrixD();
    Terr = new TMatrixD();
    SVD = new TDecompSVD();
    rand = new TRandom3(0);
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
    delete Terr;
    delete SVD;
    delete rand;
}

void Pidrix::SetTarget(TH2 *target) {
    maxm = target->GetNbinsY();
    maxn = target->GetNbinsX();

    T->ResizeTo(maxm, maxn);
    Terr->ResizeTo(maxm, maxn);

    integral = 0;
    for(int i = 0; i < maxm; i++) {
        for(int j = 0; j < maxn; j++) {
            (*T)[i][j] = target->GetBinContent(i+1, j+1);
            (*Terr)[i][j] = target->GetBinError(i+1, j+1);
            integral += (*T)[i][j];
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

const TVectorD& Pidrix::SVDSigma() {
    return SVD->GetSig();
}

double Pidrix::RandomUniform(double min, double max) {
    return rand->Uniform(min, max);
}
