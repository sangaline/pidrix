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

};

#endif
