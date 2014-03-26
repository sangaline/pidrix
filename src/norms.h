#ifndef norms_h
#define norms_h

#include "TMatrixDfwd.h"
#include "TVectorDfwd.h"
class Pidrix;
class TH1;

namespace Norms {
    double Euclidian(Pidrix *P);
    double Euclidian(const TMatrixD* A, const TMatrixD* B);
    double Euclidian(const TVectorD* A, const TVectorD* B);
    double Euclidian(const TH1* A, const TH1* B);

    //note that these are only truly KBD when both A and B are normalized
    double KullbackLeibler(Pidrix *P);
    double KullbackLeibler(const TMatrixD* A, const TMatrixD* B);
    double KullbackLeibler(const TVectorD* A, const TVectorD* B);
    double KullbackLeibler(const TH1* A, const TH1* B);

    double SymmetrizedKullbackLeibler(Pidrix *P);
    double SymmetrizedKullbackLeibler(const TMatrixD* A, const TMatrixD* B);
    double SymmetrizedKullbackLeibler(const TVectorD* A, const TVectorD* B);
    double SymmetrizedKullbackLeibler(const TH1* A, const TH1* B);
};

#endif
