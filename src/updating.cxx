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

void Updating::MultiplicativeKL(Pidrix *P, const unsigned int iterations) {
    TMatrixD U = P->GetU();
    TMatrixD V = P->GetV();
    const TMatrixD T = P->GetT();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    for(unsigned int iteration = 0; iteration < iterations; iteration++) {
        TMatrixD A(U, TMatrixD::kMult,V);
        TMatrixD oldU(U);

        //U-Update
        for(unsigned int i = 0; i < m; i++) {
            for(unsigned int j = 0; j < rank; j++) {
                double sum = 0;
                for(unsigned int mu = 0; mu < n; mu++) {
                    if( A[i][mu] == 0) { continue; } //avoid division by zero
                    sum += V[j][mu]*T[i][mu]/A[i][mu];
                }
                U[i][j] *= sum;

                sum = 0;
                for(unsigned int mu = 0; mu < n; mu++) {
                    sum += V[j][mu];
                }
                if(sum == 0) {
                    U[i][j] = 0; 
                }
                else {
                    U[i][j] /= sum;
                }
            }
        }

        //V-Update
        for(unsigned int i = 0; i < rank; i++) {
            for(unsigned int j = 0; j < n; j++) {
                double sum = 0;
                for(unsigned int mu = 0; mu < m; mu++) {
                    if( A[mu][j] == 0) { continue; } //avoid division by zero
                    sum += oldU[mu][i]*T[mu][j]/A[mu][j];
                }
                V[i][j] *= sum;

                sum = 0;
                for(unsigned int mu = 0; mu < m; mu++) {
                    sum += oldU[mu][i];
                }
                if(sum == 0) {
                    V[i][j] = 0; 
                }
                else {
                    V[i][j] /= sum;
                }
            }
        }
    }

    P->SetU(U);
    P->SetV(V);
}
