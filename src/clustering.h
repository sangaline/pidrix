#ifndef clustering_h
#define clustering_h

#include "TMatrixDfwd.h"
class Pidrixter;
class Pidrix;
namespace Norms { double SymmetrizedKullbackLeibler(const TMatrixD* A, const TMatrixD* B); }

namespace Clustering {
    void KMeans(Pidrixter *PXT, unsigned int iterations = 2, double (*norm)(const TMatrixD*, const TMatrixD*) = Norms::SymmetrizedKullbackLeibler);

    void Reorder(Pidrix *P, unsigned int* order);
    void Reorder(Pidrixter *PXT, unsigned int* order);
};

#endif
