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

double Norms::MeansAndYields(const TMatrixD* T, const TMatrixD* A) {
    double sumA_x = 0, sumT_x = 0, sum2T_x = 0;
    double sumA_y = 0, sumT_y = 0, sum2T_y = 0;
    double yieldT = 0, yieldA = 0;
    const unsigned int m = T->GetNrows();
    const unsigned int n = T->GetNcols();
    for(unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0; j < n; j++) {
            const double t = (*T)[i][j];
            const double a = (*A)[i][j];

            sumA_x += a*double(i);
            sumA_y += a*double(j);
            yieldA += a;

            sumT_x += t*double(i);
            sum2T_x += t*double(i*i);
            sumT_y += t*double(j);
            sum2T_y += t*double(j*j);
            yieldT += t;
        }
    }
    double meanA_x = sumA_x/yieldA;
    double meanA_y = sumA_y/yieldA;

    double meanT_x = sumT_x/yieldT;
    double meanT_y = sumT_y/yieldT;
    double varianceT_x = (sum2T_x/yieldT) - meanT_x*meanT_x;
    double varianceT_y = (sum2T_y/yieldT) - meanT_y*meanT_y;

    double chi2 = 0;
    chi2 += pow( meanT_x - meanA_x, 2)/varianceT_x;
    chi2 += pow( meanT_y - meanA_y, 2)/varianceT_y;
    chi2 += pow( yieldT - yieldA, 2)/yieldT;
    chi2 /= 3.0;

    return chi2;
}

double Norms::SymmetrizedMeansAndYields(const TMatrixD* A, const TMatrixD* B) {
    return 0.5*(MeansAndYields(A, B) + MeansAndYields(B, A));
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
