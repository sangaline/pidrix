#include "plotting.h"

#include "Pidrix.h"

#include "TGraph.h"
#include "TH2D.h"

#include "TVectorD.h"
#include "TMatrixD.h"

using namespace Plotting;

TGraph* Plotting::SVGraph(Pidrix *P, TGraph *t) {
    if(t == 0) {
        t = new TGraph();
    }
    else {
        t->Set(0);
    }

    const TVectorD& S = P->SVDSigma();
    const int Slen = S.GetNoElements();

    for(int i = 0; i < Slen; i++) {
        t->SetPoint(i, i+1, S[i]);
    }

    return t;
}

TH2D* Plotting::Approximation(Pidrix* P, TH2D* h, const char* name) {
    if(h == 0) {
        h = new TH2D(name, "Current Approximation;x;y", 
            P->Columns(), P->LowX(), P->HighX(),
            P->Rows(), P->LowY(), P->HighY());
    }
    else {
        h->Reset();
    }
    const TMatrixD& U = P->GetU();
    const TMatrixD& V = P->GetV();
    TMatrixD A(U,TMatrixD::kMult,V);

    const int m = P->Rows();
    const int n = P->Columns();
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            h->SetBinContent(i+1, j+1, A[i][j]);
        }
    }

    return h;
}

TH2D* Plotting::Target(Pidrix* P, TH2D* h, const char* name) {
    if(h == 0) {
        h = new TH2D(name, "Original Distribution;x;y", 
            P->Columns(), P->LowX(), P->HighX(),
            P->Rows(), P->LowY(), P->HighY());
    }
    else {
        h->Reset();
    }
    const TMatrixD& T = P->GetT();
    const TMatrixD& E = P->GetE();

    const int m = P->Rows();
    const int n = P->Columns();
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            h->SetBinContent(i+1, j+1, T[i][j]);
            h->SetBinError(i+1, j+1, T[i][j]);
        }
    }

    return h;
}
