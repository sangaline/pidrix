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
