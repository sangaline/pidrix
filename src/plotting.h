#ifndef plotting_h
#define plotting_h

#include "TMatrixDfwd.h"
class Pidrixter;
class Pidrix;
class TGraph;
class TH2D;
class TH1D;

namespace Norms { double SymmetrizedKullbackLeibler(const TMatrixD* A, const TMatrixD* B); }

namespace Plotting {
    TGraph* SVGraph(const Pidrix* P, TGraph *t = 0);
    TH2D* Approximation(const Pidrix* P, TH2D* h = 0);
    TH2D* Target(const Pidrix* P, TH2D* h = 0);
    TH1D* DistributionX(const Pidrix* P, unsigned int vector, TH1D* h = 0);
    TH1D* DistributionY(const Pidrix* P, unsigned int vector, TH1D* h = 0);
    TH2D* DistributionXY(const Pidrix* P, unsigned int vector, TH2D* h = 0);

    //versions that take in a Pidrixter and compute the errors
    TH2D* DistributionXY(Pidrixter* PDX, unsigned int vector, TH2D* h = 0);
    TH1D* DistributionX(Pidrixter* PDX, unsigned int vector, TH1D* h = 0);
    TH1D* DistributionY(Pidrixter* PDX, unsigned int vector, TH1D* h = 0);

    //only for rank 2 Pidrixes
    TGraph** Clusters(Pidrixter* PXT, TGraph** t = 0, double (*norm)(const TMatrixD*, const TMatrixD*) = Norms::SymmetrizedKullbackLeibler);

};

#endif
