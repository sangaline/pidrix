#include "plotting.h"

#include "Pidrix.h"

#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

#include "TVectorD.h"
#include "TMatrixD.h"

using namespace Plotting;

TGraph* Plotting::SVGraph(const Pidrix *P, TGraph *t) {
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

TH2D* Plotting::Approximation(const Pidrix* P, TH2D* h, const char* name) {
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
            h->SetBinContent(j+1, i+1, A[i][j]);
        }
    }

    return h;
}

TH2D* Plotting::Target(const Pidrix* P, TH2D* h, const char* name) {
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

    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    for(unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0; j < n; j++) {
            h->SetBinContent(j+1, i+1, T[i][j]);
            h->SetBinError(j+1, i+1, T[i][j]);
        }
    }

    return h;
}

TH1D* Plotting::DistributionX(const Pidrix* P, unsigned int vector, TH1D* h) {
    if(h == 0) {
        h = new TH1D("", "Distribution X;x", 
            P->Columns(), P->LowX(), P->HighX());
    }
    else {
        h->Reset();
    }

    const TMatrixD& U = P->GetU();
    const TMatrixD& V = P->GetV();

    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();

    double U_contribution = 0;
    for(unsigned int i = 0; i < m; i++) {
        U_contribution += U[i][vector];
    }
    for(unsigned int j = 0; j < n; j++) {
        h->SetBinContent(j+1, V[vector][j]*U_contribution);
    }

    return h;
}

TH1D* Plotting::DistributionY(const Pidrix* P, unsigned int vector, TH1D* h) {
    if(h == 0) {
        h = new TH1D("", "Distribution Y;y", 
            P->Rows(), P->LowY(), P->HighY());
    }
    else {
        h->Reset();
    }

    const TMatrixD& U = P->GetU();
    const TMatrixD& V = P->GetV();

    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();

    double V_contribution = 0;
    for(unsigned int j = 0; j < n; j++) {
        V_contribution += V[vector][j];
    }
    for(unsigned int i = 0; i < m; i++) {
        h->SetBinContent(i+1, U[i][vector]*V_contribution);
    }

    return h;
}
