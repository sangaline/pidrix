#include "plotting.h"

#include "Pidrix.h"
#include "Pidrixter.h"
#include "norms.h"

#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

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

TH2D* Plotting::Approximation(const Pidrix* P, TH2D* h) {
    if(h == 0) {
        h = new TH2D("", "Current Approximation;x;y", 
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
            h->SetBinError(j+1, i+1, TMath::Sqrt(A[i][j]));
        }
    }

    return h;
}

TH2D* Plotting::Target(const Pidrix* P, TH2D* h) {
    if(h == 0) {
        h = new TH2D("", "Original Distribution;x;y", 
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
            h->SetBinError(j+1, i+1, E[i][j]);
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
        h->SetBinError(j+1, TMath::Sqrt(V[vector][j]*U_contribution));
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
        h->SetBinError(i+1, TMath::Sqrt(U[i][vector]*V_contribution));
    }

    return h;
}

TH2D* Plotting::DistributionXY(const Pidrix* P, unsigned int vector, TH2D* h) {
    if(h == 0) {
        h = new TH2D("", "Distribution XY;x;y", 
            P->Columns(), P->LowX(), P->HighX(),
            P->Rows(), P->LowY(), P->HighY());
    }
    else {
        h->Reset();
    }
    const TMatrixD& U = P->GetU();
    const TMatrixD& V = P->GetV();

    const int m = P->Rows();
    const int n = P->Columns();
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            h->SetBinContent(j+1, i+1, U[i][vector]*V[vector][j]);
            h->SetBinError(j+1, i+1, TMath::Sqrt(U[i][vector]*V[vector][j]));
        }
    }

    return h;
}

TGraph** Plotting::Clusters(Pidrixter* PXT, TGraph** t, double (*norm)(const TMatrixD*, const TMatrixD*)) {
    Pidrix *P = PXT->Member(0);
    const unsigned int rank = P->Rank();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    TMatrixD Ucolumn(m,1), Vrow(1,n);
    TMatrixD Means[2] = {TMatrixD(m,n), TMatrixD(m,n)};
    Means[0].Zero();
    Means[1].Zero();

    unsigned int count = 0;
    for(unsigned int p = 0; p <  PXT->Members(); p++) {
        P = PXT->Member(p);
        TMatrixD U = P->GetU();
        TMatrixD V = P->GetV();
        for(int mu = 0; mu < 2; mu++) {
            TMatrixDColumn(Ucolumn, 0) = TMatrixDColumn(U, mu);
            TMatrixDRow(Vrow, 0) = TMatrixDRow(V, mu);
            Means[mu] += Ucolumn*Vrow;
        }
        count++;
    }
    Means[0] *= 1.0/double(count);
    Means[1] *= 1.0/double(count);

    if(t == 0) {
        t = new TGraph* [3];
        t[0] = new TGraph();
        t[1] = new TGraph();
        t[2] = new TGraph();
    }
    else {
        t[0]->Set(0);
        t[1]->Set(0);
        t[2]->Set(0);
    }

    double d = norm(&Means[0], &Means[1]);
    TMatrixD UV(m,n);
    for(unsigned int p = 0; p <  PXT->Members(); p++) {
        P = PXT->Member(p);
        TMatrixD U = P->GetU();
        TMatrixD V = P->GetV();
        //mu is the class
        for(int mu = 0; mu < 2; mu++) {
            TMatrixDColumn(Ucolumn, 0) = TMatrixDColumn(U, mu);
            TMatrixDRow(Vrow, 0) = TMatrixDRow(V, mu);
            UV = Ucolumn*Vrow;
            double r0 = norm(&UV, &Means[0]);
            double r1 = norm(&UV, &Means[1]);
            t[mu]->SetPoint(p, r0, r1);
            t[2]->SetPoint(p*2+mu, r0, r1);
        }
    }
    return t;
}
