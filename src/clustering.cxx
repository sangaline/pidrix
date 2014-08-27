#include "clustering.h"

#include "Pidrixter.h"
#include "Pidrix.h"
#include "norms.h"

#include "TMatrixD.h"
#include "TMath.h"

using namespace Clustering;

void Clustering::Reorder(Pidrix *P, unsigned int* order) {
    const unsigned int rank = P->Rank();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    TMatrixD oldU = P->GetU();
    TMatrixD oldV = P->GetV();
    TMatrixD U(m, rank);
    TMatrixD V(rank, n);

    for(unsigned int mu = 0; mu < rank; mu++) {
        TMatrixDColumn(U, mu) = TMatrixDColumn(oldU, order[mu]);
        TMatrixDRow(V, mu) = TMatrixDRow(oldV, order[mu]);
    }
    P->SetU(U);
    P->SetV(V);
}

void Clustering::Reorder(Pidrixter *PXT, unsigned int* order) {
    for(unsigned int p = 0; p <  PXT->Members(); p++) {
        Reorder(PXT->Member(p), order);
    }
}

void Clustering::KMeans(Pidrixter *PXT, unsigned int iterations, double (*norm)(const TMatrixD*, const TMatrixD*)) {
    Pidrix *P = PXT->Member(0);
    const unsigned int rank = P->Rank();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    TMatrixD Ucolumn(m,1), Vrow(1,n);

    TMatrixD* Sum = new TMatrixD [rank];
    TMatrixD* Current = new TMatrixD [rank];
    TMatrixD* Mean = new TMatrixD [rank];
    unsigned int* order = new unsigned int[rank];
    unsigned int count = 0;
    TMatrixD Distances(rank, rank);
    for(unsigned int i = 0; i < rank; i++) {
        Sum[i].ResizeTo(m,n);
        Current[i].ResizeTo(m,n);
        Mean[i].ResizeTo(m,n);

        Sum[i] = TMatrixD(m,n);
        Sum[i].Zero();
    }

    for(unsigned int iteration = 0; iteration < iterations; iteration++) {
        for(unsigned int p = 0; p <  PXT->Members(); p++) {
            P = PXT->Member(p);
            TMatrixD U = P->GetU();
            TMatrixD V = P->GetV();

            //update the means and the current matrices
            for(unsigned int mu = 0; mu < rank; mu++) {
                TMatrixDColumn(Ucolumn, 0) = TMatrixDColumn(U, mu);
                TMatrixDRow(Vrow, 0) = TMatrixDRow(V, mu);
                Current[mu] = Ucolumn*Vrow;

                if(iteration > 0 || p > 0) {
                    Mean[mu] = Sum[mu]*double(1.0/double(count));
                }
                else {
                    Sum[mu] += Current[mu];
                }
            }
            if( iteration == 0 && p == 0 ) {
                count++;
                continue;
            }

            //compute the distance matrix
            for(unsigned int i = 0; i < rank; i++) {
                for(unsigned int j = 0; j < rank; j++) {
                    Distances[i][j] = norm(&Current[i], &Mean[j]);
                }
            }
            //find the optimal matchings
            for(unsigned int mu = 0; mu < rank; mu++) {
                double minval = -1;
                unsigned int mini, minj;
                //find the closest remaining match
                for(unsigned int i = 0; i < rank; i++) {
                    for(unsigned int j = 0; j < rank; j++) {
                        //norms are positive semi-definite so we use negative values as invalid
                        if(Distances[i][j] >= 0) {
                            if(Distances[i][j] < minval || minval < 0) {
                                mini = i;
                                minj = j;
                                minval = Distances[i][j];
                            }
                        }
                    }
                }
                //mind i vs j here, we want order[0] to be the index of what the first element 
                //should be switched to in order to match the means
                order[minj] = mini;
                //invalidate matches that were just eliminated
                for(unsigned int nu = 0; nu < rank; nu++) {
                    Distances[mini][nu] = -1;
                    Distances[nu][minj] = -1;
                }
            }

            //update the sum
            count++;
            for(unsigned int mu = 0; mu < rank; mu++) {
                Sum[mu] += Current[order[mu]];
            }

            //reorder this pidrix to match the current means if this is the last iteration
            if(iteration+1 == iterations) {
                Reorder(P, order);
            }
        }
    }

    //reorder them from largest to smallest yields
    double* yields = new double[rank];
    for(unsigned int mu = 0; mu < rank; mu++) {
        yields[mu] = Mean[mu].Norm1();
    }
    TMath::Sort(rank, yields, order, true);
    Reorder(PXT, order);
    
    //cleanup
    delete [] Sum;
    delete [] Current;
    delete [] Mean;
    delete [] order;
}
