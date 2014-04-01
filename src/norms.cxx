#include "norms.h"
#include "Pidrix.h"

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH1.h"

#include "math.h"

using namespace Norms;

double Norms::Euclidian(Pidrix *P) {
    const TMatrixD A = P->GetU()*P->GetV();
    return Euclidian(&P->GetT(), &A);
}

double Norms::Euclidian(const TMatrixD* A, const TMatrixD* B) {
    return sqrt(((*A)-(*B)).E2Norm());
}

double Norms::Euclidian(const TVectorD* A, const TVectorD* B) {
    return sqrt(((*A)-(*B)).Norm2Sqr());
}

double Norms::Euclidian(const TH1* A, const TH1* B) {
    const unsigned int xbins = A->GetNbinsX();
    const unsigned int ybins = A->GetNbinsY();
    double sum2 = 0;
    for(unsigned int i = 1; i <= xbins; i++) {
        for(unsigned int j = 1; j <= ybins; j++) {
            sum2 += pow(A->GetBinContent(i,j) - B->GetBinContent(i,j), 2);
        }
    }
    return sqrt( sum2/double(xbins*ybins) );
}


double Norms::KullbackLeibler(Pidrix *P) {
    const TMatrixD A = P->GetU()*P->GetV();
    return KullbackLeibler(&P->GetT(), &A);
}

double Norms::KullbackLeibler(const TMatrixD* A, const TMatrixD* B) {
    const unsigned int m = A->GetNrows();
    const unsigned int n = A->GetNcols();
    double sum = 0;
    for(unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0; j < n; j++) {
            const double a = (*A)[i][j];
            const double b = (*B)[i][j];
            if(b > 0 && a > 0) {
                sum += a*log(a/b);
            }
            sum += - a + b;
        }
    }
    return sum;
}

double Norms::KullbackLeibler(const TVectorD* A, const TVectorD* B) {
    const unsigned int m = A->GetNoElements();
    double sum = 0;
    for(unsigned int i = 0; i < m; i++) {
        const double a = (*A)[i];
        const double b = (*B)[i];
        if(b > 0) {
            sum += a*log(a/b) - a + b;
        }
    }
    return sum;
}

double Norms::KullbackLeibler(const TH1* A, const TH1* B) {
    const unsigned int xbins = A->GetNbinsX();
    const unsigned int ybins = A->GetNbinsY();
    double sum = 0;
    for(unsigned int i = 1; i <= xbins; i++) {
        for(unsigned int j = 1; j <= ybins; j++) {
            const double a = A->GetBinContent(i,j);
            const double b = B->GetBinContent(i,j);
            if(b > 0) {
                sum += a*log(a/b) - a + b;
            }
        }
    }
    return sum;
}

double Norms::SymmetrizedKullbackLeibler(Pidrix *P) {
    const TMatrixD A = P->GetU()*P->GetV();
    return SymmetrizedKullbackLeibler(&P->GetT(), &A);
}

double Norms::SymmetrizedKullbackLeibler(const TMatrixD* A, const TMatrixD* B) {
    return 0.5*(KullbackLeibler(A, B) + KullbackLeibler(B, A));
}

double Norms::SymmetrizedKullbackLeibler(const TVectorD* A, const TVectorD* B) {
    return 0.5*(KullbackLeibler(A, B) + KullbackLeibler(B, A));
}

double Norms::SymmetrizedKullbackLeibler(const TH1* A, const TH1* B) {
    return 0.5*(KullbackLeibler(A, B) + KullbackLeibler(B, A));
}

double Norms::ChiSquared(Pidrix *P, bool per_ndf) {
    const TMatrixD& T = P->GetT();
    const TMatrixD& E = P->GetE();
    const TMatrixD& A = P->GetU()*P->GetV();
    const unsigned int m = A.GetNrows();
    const unsigned int n = A.GetNcols();
    double sum = 0, ndf = 0;
    for(unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0; j < n; j++) {
            const double t = T[i][j];
            const double e = E[i][j];
            const double a = A[i][j];
            if(t > 0 && e > 0) {
                sum += pow((t-a)/e, 2);
                ndf += 1;
            }
        }
    }
    if(per_ndf) {
        return sum/ndf;
    }
    else {
        return sum;
    }
}

double Norms::ChiSquared(const TMatrixD* T, const TMatrixD* A) {
    const unsigned int m = A->GetNrows();
    const unsigned int n = A->GetNcols();
    double sum = 0;
    for(unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0; j < n; j++) {
            if( (*T)[i][j] > 0) {
                sum += pow((*T)[i][j]-(*A)[i][j],2)/(*T)[i][j];
            }
        }
    }
    return sum;
}

double Norms::SymmetrizedChiSquared(const TMatrixD* A, const TMatrixD* B) {
    return 0.5*(ChiSquared(A,B)+ChiSquared(B,A));
}

double Norms::Yields(const TMatrixD* T, const TMatrixD* A) {
    return fabs(T->Sum()-A->Sum())/sqrt(T->Sum());
}

double Norms::SymmetrizedYields(const TMatrixD* A, const TMatrixD* B) {
    return fabs(B->Sum()-A->Sum())/sqrt(pow(A->Sum(),2)+pow(B->Sum(), 2));
}
