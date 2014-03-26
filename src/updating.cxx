#include "updating.h"

#include "Pidrix.h"

#include "TMatrixD.h"

using namespace Updating;

//Sets division by zero to zero rather than naan
void CleanElementDiv(TMatrixD& A, TMatrixD&B) {
    const unsigned int m = A.GetNrows();
    const unsigned int n = A.GetNcols();
    for(unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0; j < n; j++) {
            if(B[i][j] == 0) {
                A[i][j] = 0;
            }
            else {
                A[i][j] /= B[i][j];
            }
        }
    }
}

void Updating::MultiplicativeEuclidian(Pidrix *P, const unsigned int iterations) {
    TMatrixD U = P->GetU();
    TMatrixD V = P->GetV();
    const TMatrixD T = P->GetT();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    for(unsigned int iteration = 0; iteration < iterations; iteration++) {
        //element-wise: new U = U*A1/B1
        //element-wise: new V = V*A2/B2
        TMatrixD A1(T, TMatrixD::kMultTranspose,V);
        TMatrixD B1(U,TMatrixD::kMult,TMatrixD(V,TMatrixD::kMultTranspose,V));
        TMatrixD A2(U, TMatrixD::kTransposeMult,T);
        TMatrixD B2(TMatrixD(U,TMatrixD::kTransposeMult,U),TMatrixD::kMult,V);

        CleanElementDiv(A1, B1);
        ElementMult(U, A1);

        CleanElementDiv(A2, B2);
        ElementMult(V, A2);
    }
    P->SetU(U);
    P->SetV(V);
}
