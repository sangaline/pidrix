#include "quantifying.h"

#include "Pidrix.h"

#include "TVectorD.h"
#include "TMatrixD.h"

#include "math.h"

using namespace Quantifying;

double Quantifying::Yield(Pidrix *P, unsigned int vector) {
    const TMatrixD& U = P->GetU();
    const TMatrixD& V = P->GetV();
    const unsigned int maxm = P->Rows();
    const unsigned int maxn = P->Columns();
    double Usum = 0;
    for(unsigned int i = 0; i < maxm; i++) {
        Usum += U[i][vector];
    }
    double Vsum = 0;
    for(unsigned int j = 0; j < maxn; j++) {
        Vsum += V[vector][j];
    }
    return Vsum*Usum;
}

double Quantifying::MeanX(Pidrix *P, unsigned int vector) {
    const TMatrixD& V = P->GetV();
    const unsigned int maxn = P->Columns();
    double xsum = 0, Vsum = 0;
    for(unsigned int j = 0; j < maxn; j++) {
        Vsum += V[vector][j];
        xsum += V[vector][j]*X(P, j);
    }
    return xsum/Vsum;
}

double Quantifying::MeanY(Pidrix *P, unsigned int vector) {
    const TMatrixD& U = P->GetU();
    const unsigned int maxm = P->Rows();
    double ysum = 0, Usum = 0;
    for(unsigned int i = 0; i < maxm; i++) {
        Usum += U[i][vector];
        ysum += U[i][vector]*Y(P, i);
    }
    return ysum/Usum;
}

double Quantifying::VarianceX(Pidrix *P, unsigned int vector) {
    const TMatrixD& V = P->GetV();
    const unsigned int maxn = P->Columns();
    double xsum = 0, x2sum = 0, Vsum = 0;
    for(unsigned int j = 0; j < maxn; j++) {
        const double x = X(P, j);
        Vsum += V[vector][j];
        xsum += V[vector][j]*x;
        x2sum += V[vector][j]*x*x;
    }
    return (x2sum/Vsum) - ((xsum/Vsum)*(xsum/Vsum));
}

double Quantifying::VarianceY(Pidrix *P, unsigned int vector) {
    const TMatrixD& U = P->GetU();
    const unsigned int maxm = P->Rows();
    double ysum = 0, ysum2 = 0, Usum = 0;
    for(unsigned int i = 0; i < maxm; i++) {
        const double y = Y(P, i);
        Usum += U[i][vector];
        ysum += U[i][vector]*y;
        ysum2 += U[i][vector]*y*y;
    }
    return (ysum2/Usum) - ((ysum/Usum)*(ysum/Usum));
}

double Quantifying::StandardDeviationX(Pidrix *P, unsigned int vector) {
    return sqrt(VarianceX(P, vector));
}

double Quantifying::StandardDeviationY(Pidrix *P, unsigned int vector) {
    return sqrt(VarianceY(P, vector));
}

double Quantifying::X(Pidrix *P, unsigned int j) {
    const double xlow = P->LowX();
    const double xhigh = P->HighX();
    const double delta = (xhigh - xlow)/double(P->Columns());
    return xlow + (double(j) + 0.5)*delta;
}

double Quantifying::Y(Pidrix *P, unsigned int i) {
    const double ylow = P->LowY();
    const double yhigh = P->HighY();
    const double delta = (yhigh - ylow)/double(P->Rows());
    return ylow + (double(i) + 0.5)*delta;
}
