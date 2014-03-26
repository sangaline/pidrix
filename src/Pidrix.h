#ifndef Pidrix_h
#define Pidrix_h

#include "TNamed.h"

#include "TMatrixDfwd.h"
#include "TVectorDfwd.h"
class TDecompSVD;

class TH2;
class TGraph;
class TRandom;

class Pidrix : public TNamed {
    TMatrixD *T, *Terr; // T for target
    double integral, xlow, xhigh, ylow, yhigh;
    TMatrixD *U, *V; // UV ~= T
    unsigned int maxm, maxn;
    unsigned int rank;

    //Various things that helper functions will rely on
    TDecompSVD *SVD;
    TRandom *rand;

  public:
    Pidrix();
    ~Pidrix();

    void SetTarget(TH2 *target);

    TGraph* SVGraph(TGraph *t = 0);

    unsigned int Rank() { return rank; }
    double Integral() { return integral; }
    unsigned int Rows() { return maxm; }
    unsigned int Columns() { return maxn; }
    double LowX() { return xlow; }
    double HighX() { return xhigh; }
    double LowY() { return ylow; }
    double HighY() { return yhigh; }

    unsigned int SetRank(unsigned int new_rank);
    double SetU(unsigned i, unsigned j, double value);
    double SetV(unsigned i, unsigned j, double value);
    double GetU(unsigned i, unsigned j);
    double GetV(unsigned i, unsigned j);
    const TMatrixD& GetU() { return (*U); }
    const TMatrixD& GetV() { return (*V); }

    double RandomUniform(double min, double max);

    const TVectorD& SVDSigma();

    ClassDef(Pidrix,1) //PID Matrix Factorization
};

#endif
