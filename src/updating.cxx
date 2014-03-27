#include "updating.h"

#include "Pidrix.h"
#include "quantifying.h"

#include "TMatrixD.h"

#include "math.h"

using namespace Updating;

//Dampens division by zero
void CleanElementDiv(TMatrixD& A, TMatrixD&B, double epsilon = 1e-16) {
    const unsigned int m = A.GetNrows();
    const unsigned int n = A.GetNcols();
    for(unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0; j < n; j++) {
            A[i][j] /= B[i][j] + epsilon;
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
        //element-wise: new U = U*A/B
        TMatrixD A(T, TMatrixD::kMultTranspose,V);
        TMatrixD B(U,TMatrixD::kMult,TMatrixD(V,TMatrixD::kMultTranspose,V));

        CleanElementDiv(A, B);
        ElementMult(U, A);

        //element-wise: new V = V*C/D
        TMatrixD C(U, TMatrixD::kTransposeMult,T);
        TMatrixD D(TMatrixD(U,TMatrixD::kTransposeMult,U),TMatrixD::kMult,V);
        CleanElementDiv(C, D);
        ElementMult(V, C);
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

        A = U*V;
        //V-Update
        for(unsigned int i = 0; i < rank; i++) {
            for(unsigned int j = 0; j < n; j++) {
                double sum = 0;
                for(unsigned int mu = 0; mu < m; mu++) {
                    if( A[mu][j] == 0) { continue; } //avoid division by zero
                    sum += U[mu][i]*T[mu][j]/A[mu][j];
                }
                V[i][j] *= sum;

                sum = 0;
                for(unsigned int mu = 0; mu < m; mu++) {
                    sum += U[mu][i];
                }
                V[i][j] /= sum;
            }
        }
    }

    P->SetU(U);
    P->SetV(V);
}

void Updating::Normalize(Pidrix *P) {
    TMatrixD U = P->GetU();
    TMatrixD V = P->GetV();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    const double integral = P->Integral();
    const double scale = P->Integral() / (U*V).Norm1();
    for(int vector = 0; vector < rank; vector++) {
        double Usum = 0;
        for(unsigned int i = 0; i < m; i++) {
            Usum += U[i][vector];
        }
        double Vsum = 0;
        for(unsigned int j = 0; j < n; j++) {
            Vsum += V[vector][j];
        }
        const double Uscale = sqrt(scale*Vsum/Usum);
        const double Vscale = sqrt(scale*Usum/Vsum);

        for(unsigned int i = 0; i < m; i++) {
            U[i][vector] *= Uscale;
        }
        for(unsigned int j = 0; j < n; j++) {
            V[vector][j] *= Vscale;
        }
    }

    P->SetU(U);
    P->SetV(V);
}

void Updating::ScaleY(Pidrix *P, double factor) {
    const TMatrixD& oldU = P->GetU();
    TMatrixD U = oldU;

    const unsigned int m = P->Rows();
    const unsigned int rank = P->Rank();

    const double ylow = P->LowY();
    const double yhigh = P->HighY();
    const double ydelta = (yhigh - ylow)/double(P->Rows());

    for(int vector = 0; vector < rank; vector++) {
        double meany = Quantifying::MeanY(P, vector);
        for(unsigned int i = 0; i < m; i++) {
            //after transformation
            double lowedge = ylow + double(i)*ydelta;
            double highedge = lowedge + ydelta;
            //before transformation
            lowedge = meany + ((lowedge - meany)/factor);
            highedge = meany + ((highedge - meany)/factor);
            //which bins would they be in?
            const int lowbin = floor((lowedge - ylow)/ydelta);
            const int highbin = floor((highedge - ylow)/ydelta);
            double value = 0;
            if(lowbin == highbin) {
                if(lowbin >= 0 && lowbin < m) {
                    const double bin_fraction = (highedge-lowedge)/ydelta;
                    value = bin_fraction*oldU[lowbin][vector];
                }
            }
            else {
                if(lowbin >= 0 && lowbin < m) {
                    const double lowbin_highedge = ylow + double(lowbin+1)*ydelta;
                    const double bin_fraction = (lowbin_highedge-lowedge)/ydelta;
                    value += bin_fraction*oldU[lowbin][vector];
                }
                if(highbin >= 0 && highbin < m) {
                    const double highbin_lowedge = ylow + double(highbin)*ydelta;
                    const double bin_fraction = (highedge - highbin_lowedge)/ydelta;
                    value += bin_fraction*oldU[highbin][vector];
                }
                for(int bin = lowbin+1; bin < highbin; bin++) {
                    if(bin >= 0 && bin < m) {
                        value += oldU[bin][vector];
                    }
                }
            }
            U[i][vector] = value;
        }
    }
    P->SetU(U);
}

void Updating::ScaleX(Pidrix *P, double factor) {
    const TMatrixD& oldV = P->GetV();
    TMatrixD V = oldV;

    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    const double xlow = P->LowX();
    const double xhigh = P->HighX();
    const double xdelta = (xhigh - xlow)/double(P->Columns());

    for(int vector = 0; vector < rank; vector++) {
        double meanx = Quantifying::MeanX(P, vector);
        for(unsigned int j = 0; j < n; j++) {
            //after transformation
            double lowedge = xlow + double(j)*xdelta;
            double highedge = lowedge + xdelta;
            //before transformation
            lowedge = meanx + ((lowedge - meanx)/factor);
            highedge = meanx + ((highedge - meanx)/factor);
            //which bins would they be in?
            const int lowbin = floor((lowedge - xlow)/xdelta);
            const int highbin = floor((highedge - xlow)/xdelta);
            double value = 0;
            if(lowbin == highbin) {
                if(lowbin >= 0 && lowbin < n) {
                    const double bin_fraction = (highedge-lowedge)/xdelta;
                    value = bin_fraction*oldV[vector][lowbin];
                }
            }
            else {
                if(lowbin >= 0 && lowbin < n) {
                    const double lowbin_highedge = xlow + double(lowbin+1)*xdelta;
                    const double bin_fraction = (lowbin_highedge-lowedge)/xdelta;
                    value += bin_fraction*oldV[vector][lowbin];
                }
                if(highbin >= 0 && highbin < n) {
                    const double highbin_lowedge = xlow + double(highbin)*xdelta;
                    const double bin_fraction = (highedge - highbin_lowedge)/xdelta;
                    value += bin_fraction*oldV[vector][highbin];
                }
                for(int bin = lowbin+1; bin < highbin; bin++) {
                    if(bin >= 0 && bin < n) {
                        value += oldV[vector][bin];
                    }
                }
            }
            V[vector][j] = value;
        }
    }
    P->SetV(V);
}

void Updating::Scale(Pidrix *P, double factor) {
    ScaleX(P, factor);
    ScaleY(P, factor);
}

void Updating::AddNoiseX(Pidrix *P, double fraction) {
    const TMatrixD& oldV = P->GetV();
    TMatrixD V = oldV;

    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    for(int vector = 0; vector < rank; vector++) {
        double Vsum = 0;
        for(unsigned int j = 0; j < n; j++) {
            Vsum += oldV[vector][j];
        }
        double max_noise = Vsum*fraction*2.0/double(n);
        for(unsigned int j = 0; j < n; j++) {
            V[vector][j] += P->RandomUniform(0, max_noise);
        }
    }
    P->SetV(V);
}

void Updating::AddNoiseY(Pidrix *P, double fraction) {
    const TMatrixD& oldU = P->GetU();
    TMatrixD U = oldU;

    const unsigned int m = P->Rows();
    const unsigned int rank = P->Rank();

    for(int vector = 0; vector < rank; vector++) {
        double Usum = 0;
        for(unsigned int i = 0; i < m; i++) {
            Usum += oldU[i][vector];
        }
        double max_noise = Usum*fraction*2.0/double(m);
        for(unsigned int i = 0; i < m; i++) {
            U[i][vector] += P->RandomUniform(0, max_noise);
        }
    }
    P->SetU(U);
}

void Updating::AddNoise(Pidrix *P, double fraction) {
    AddNoiseX(P, fraction);
    AddNoiseY(P, fraction);
}

void Updating::SmearX(Pidrix *P, unsigned int iterations, double neighbor_fraction) {
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    for(unsigned int iteration = 0; iteration < iterations; iteration++) {
        const TMatrixD& oldV = P->GetV();
        TMatrixD V = oldV;

        for(int vector = 0; vector < rank; vector++) {
            //edge cases
            V[vector][0] += neighbor_fraction*oldV[vector][1];
            V[vector][0] /= 1.0 + neighbor_fraction;
            V[vector][n-1] += neighbor_fraction*oldV[vector][n-2];
            V[vector][n-1] /= 1.0 + neighbor_fraction;
            for(unsigned int j = 1; j < n-1; j++) {
                V[vector][j] += neighbor_fraction*(oldV[vector][j-1]+oldV[vector][j+1]);
                V[vector][j] /= 1.0 + 2.0*neighbor_fraction;
            }
        }
        P->SetV(V);
    }
}

void Updating::SmearY(Pidrix *P, unsigned int iterations, double neighbor_fraction) {
    const unsigned int m = P->Rows();
    const unsigned int rank = P->Rank();

    for(unsigned int iteration = 0; iteration < iterations; iteration++) {
        const TMatrixD& oldU = P->GetU();
        TMatrixD U = oldU;

        for(int vector = 0; vector < rank; vector++) {
            //edge cases
            U[0][vector] += neighbor_fraction*oldU[1][vector];
            U[0][vector] /= 1.0 + neighbor_fraction;
            U[m-1][vector] += neighbor_fraction*oldU[m-2][vector];
            U[m-1][vector] /= 1.0 + neighbor_fraction;
            for(unsigned int i = 1; i < m-1; i++) {
                U[i][vector] += neighbor_fraction*(oldU[i-1][vector]+oldU[i+1][vector]);
                U[i][vector] /= 1.0 + 2.0*neighbor_fraction;
            }
        }
        P->SetU(U);
    }
}

void Updating::Smear(Pidrix *P, unsigned int iterations, double neighbor_fraction) {
    SmearX(P, iterations, neighbor_fraction);
    SmearY(P, iterations, neighbor_fraction);
}

void Updating::ForceGaussianX(Pidrix *P) {
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    const TMatrixD& oldV = P->GetV();
    TMatrixD V = oldV;

    for(int vector = 0; vector < rank; vector++) {
        double mean = Quantifying::MeanX(P, vector);
        double variance = Quantifying::VarianceX(P, vector);
        double delta = Quantifying::X(P, 1) - Quantifying::X(P, 0);

        double Vsum = 0;
        for(unsigned int j = 0; j < n; j++) {
            Vsum += oldV[vector][j];
        }

        const double factor = Vsum/sqrt(variance*6.28318530718);
        for(unsigned int j = 0; j < n; j++) {
            V[vector][j] = factor*exp(-0.5*pow( Quantifying::X(P,j) - mean, 2)/variance);
        }
    }
    P->SetV(V);
}

void Updating::ForceGaussianY(Pidrix *P) {
    const unsigned int m = P->Rows();
    const unsigned int rank = P->Rank();

    const TMatrixD& oldU = P->GetU();
    TMatrixD U = oldU;

    for(int vector = 0; vector < rank; vector++) {
        double mean = Quantifying::MeanY(P, vector);
        double variance = Quantifying::VarianceY(P, vector);
        double delta = Quantifying::Y(P, 1) - Quantifying::Y(P, 0);

        double Usum = 0;
        for(unsigned int i = 0; i < m; i++) {
            Usum += oldU[i][vector];
        }

        const double factor = Usum/sqrt(variance*6.28318530718);
        for(unsigned int i = 0; i < m; i++) {
            U[i][vector] = factor*exp(-0.5*pow( Quantifying::Y(P,i) - mean, 2)/variance);
        }
    }
    P->SetU(U);
}

void Updating::ForceGaussian(Pidrix *P) {
    ForceGaussianY(P);
    ForceGaussianX(P);
}
