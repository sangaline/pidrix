#ifndef Pidrix_h
#define Pidrix_h

#include "TNamed.h"

#include "TMatrixDfwd.h"
#include "TVectorDfwd.h"
class TDecompSVD;

class TH2;
class TGraph;

class Pidrix : public TNamed {
    TMatrixD *T, *E; // T for target, E for Error
    double integral, xlow, xhigh, ylow, yhigh;
    TMatrixD *U, *V; // UV ~= T
    unsigned int maxm, maxn;
    unsigned int rank;

    //Various things that helper functions will rely on
    TDecompSVD *SVD;

  public:
    Pidrix();
    ~Pidrix();

    void SetTarget(TH2 *target, bool randomize = false);
    void SetTarget(const Pidrix *P, bool randomize = false);

    TGraph* SVGraph(TGraph *t = 0);

    unsigned int Rank() const { return rank; }
    double Integral() const { return integral; }
    unsigned int Rows() const { return maxm; }
    unsigned int Columns() const { return maxn; }
    double LowX() const { return xlow; }
    double HighX() const { return xhigh; }
    double LowY() const { return ylow; }
    double HighY() const { return yhigh; }

    unsigned int SetRank(unsigned int new_rank);
    double SetU(unsigned i, unsigned j, double value);
    double SetV(unsigned i, unsigned j, double value);
    double GetU(unsigned i, unsigned j);
    double GetV(unsigned i, unsigned j);
    const TMatrixD& GetU() const { return (*U); }
    const TMatrixD& GetV() const { return (*V); }
    const TMatrixD& GetT() const { return (*T); }
    const TMatrixD& GetE() const { return (*E); }
    void SetU(const TMatrixD& newU);
    void SetV(const TMatrixD& newV);

    double RandomUniform(double min, double max);

    const TVectorD& SVDSigma() const;

    ClassDef(Pidrix,1) //PID Matrix Factorization
};

#endif
