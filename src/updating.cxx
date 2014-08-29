#include "updating.h"

#include "Pidrix.h"
#include "quantifying.h"

#include "TMatrixD.h"
#include "TRandom3.h"
#include "TMath.h"

#include "math.h"

using namespace Updating;

//Dampens division by zero
void CleanElementDiv(TMatrixD& A, TMatrixD&B, double epsilon = 1e-16) {
    const unsigned int m = A.GetNrows();
    const unsigned int n = A.GetNcols();
    for(unsigned int i = 0; i < m; i++) {
        for(unsigned int j = 0; j < n; j++) {
            if(B[i][j] > 0 || epsilon > 0) {
                A[i][j] /= B[i][j] + epsilon;
            }
            else {
                A[i][j] = 0;
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

void Updating::MultiplicativeKL(Pidrix *P, const unsigned int iterations, double epsilon) {
    TMatrixD U = P->GetU();
    TMatrixD V = P->GetV();
    const TMatrixD T = P->GetT();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    TMatrixD A(m, n);
    for(unsigned int iteration = 0; iteration < iterations; iteration++) {
        A.Mult(U, V);

        //U-Update
        for(unsigned int i = 0; i < m; i++) {
            for(unsigned int j = 0; j < rank; j++) {
                double sum = 0;
                for(unsigned int mu = 0; mu < n; mu++) {
                    if(A[i][mu] > 0 || epsilon > 0) {
                        sum += V[j][mu]*T[i][mu]/(A[i][mu]+epsilon);
                    }
                }
                U[i][j] *= sum;

                sum = 0;
                for(unsigned int mu = 0; mu < n; mu++) {
                    sum += V[j][mu];
                }
                if(sum > 0 || epsilon > 0) {
                    U[i][j] /= (sum+epsilon);
                }
                else {
                    U[i][j] = 0;
                }
            }
        }

        A.Mult(U, V);
        //V-Update
        for(unsigned int i = 0; i < rank; i++) {
            for(unsigned int j = 0; j < n; j++) {
                double sum = 0;
                for(unsigned int mu = 0; mu < m; mu++) {
                    if(A[mu][j] > 0 || epsilon > 0) {
                        sum += U[mu][i]*T[mu][j]/(A[mu][j]+epsilon);
                    }
                }
                V[i][j] *= sum;

                sum = 0;
                for(unsigned int mu = 0; mu < m; mu++) {
                    sum += U[mu][i];
                }
                if(sum > 0 || epsilon > 0) {
                    V[i][j] /= (sum+epsilon);
                }
                else {
                    V[i][j] = 0;
                }
            }
        }
    }

    P->SetU(U);
    P->SetV(V);
}

void Updating::MultiplicativeKLY(Pidrix *P, const unsigned int iterations, double epsilon) {
    TMatrixD U = P->GetU();
    const TMatrixD& V = P->GetV();
    const TMatrixD& T = P->GetT();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    TMatrixD A(m, n);
    for(unsigned int iteration = 0; iteration < iterations; iteration++) {
        A.Mult(U, V);

        //U-Update
        for(unsigned int i = 0; i < m; i++) {
            for(unsigned int j = 0; j < rank; j++) {
                double sum = 0;
                for(unsigned int mu = 0; mu < n; mu++) {
                    if(A[i][mu] > 0 || epsilon > 0) {
                        sum += V[j][mu]*T[i][mu]/(A[i][mu]+epsilon);
                    }
                }
                U[i][j] *= sum;

                sum = 0;
                for(unsigned int mu = 0; mu < n; mu++) {
                    sum += V[j][mu];
                }
                if(sum > 0 || epsilon > 0) {
                    U[i][j] /= (sum+epsilon);
                }
                else {
                    U[i][j] = 0; 
                }
            }
        }
    }
    P->SetU(U);
}

void Updating::MultiplicativeKLX(Pidrix *P, const unsigned int iterations, double epsilon) {
    const TMatrixD& U = P->GetU();
    TMatrixD V = P->GetV();
    const TMatrixD& T = P->GetT();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    TMatrixD A(m, n);
    for(unsigned int iteration = 0; iteration < iterations; iteration++) {
        A.Mult(U, V);

        //V-Update
        for(unsigned int i = 0; i < rank; i++) {
            for(unsigned int j = 0; j < n; j++) {
                double sum = 0;
                for(unsigned int mu = 0; mu < m; mu++) {
                    if(A[mu][j] > 0 || epsilon > 0) {
                        sum += U[mu][i]*T[mu][j]/(A[mu][j]+epsilon);
                    }
                }
                V[i][j] *= sum;

                sum = 0;
                for(unsigned int mu = 0; mu < m; mu++) {
                    sum += U[mu][i];
                }
                if(sum > 0 || epsilon > 0) {
                    V[i][j] /= (sum+epsilon);
                }
                else {
                    V[i][j] = 0;
                }
            }
        }
    }
    P->SetV(V);
}

void Updating::Normalize(Pidrix *P) {
    TMatrixD U = P->GetU();
    TMatrixD V = P->GetV();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    const double integral = P->Integral();
    const double scale = integral / (U*V).Sum();

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
            const int lowbin = (int) floor((lowedge - ylow)/ydelta);
            const int highbin = (int) floor((highedge - ylow)/ydelta);
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
            const int lowbin = (int) floor((lowedge - xlow)/xdelta);
            const int highbin = (int)  floor((highedge - xlow)/xdelta);
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

void Updating::ForceGaussianX(Pidrix *P, unsigned int vector, double mean, double sigma) {
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    const TMatrixD& oldV = P->GetV();
    TMatrixD V = oldV;
    const double variance = sigma*sigma;
    double delta = Quantifying::X(P, 1) - Quantifying::X(P, 0);

    double Vsum = 0;
    for(unsigned int j = 0; j < n; j++) {
        Vsum += oldV[vector][j];
    }

    const double factor = Vsum/sqrt(variance*6.28318530718);
    for(unsigned int j = 0; j < n; j++) {
        V[vector][j] = factor*exp(-0.5*pow( Quantifying::X(P,j) - mean, 2)/variance);
    }
    P->SetV(V);
}

void Updating::ForceGaussianY(Pidrix *P, unsigned int vector, double mean, double sigma) {
    const unsigned int m = P->Rows();
    const unsigned int rank = P->Rank();

    const TMatrixD& oldU = P->GetU();
    TMatrixD U = oldU;
    const double variance = sigma*sigma;
    double delta = Quantifying::Y(P, 1) - Quantifying::Y(P, 0);

    double Usum = 0;
    for(unsigned int i = 0; i < m; i++) {
        Usum += oldU[i][vector];
    }

    const double factor = Usum/sqrt(variance*6.28318530718);
    for(unsigned int i = 0; i < m; i++) {
        U[i][vector] = factor*exp(-0.5*pow( Quantifying::Y(P,i) - mean, 2)/variance);
    }
    P->SetU(U);
}

void Updating::ForceGaussian(Pidrix *P) {
    ForceGaussianY(P);
    ForceGaussianX(P);
}

void Updating::ForceUnimodalX(Pidrix *P) {
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    TMatrixD V = P->GetV();

    for(int vector = 0; vector < rank; vector++) {

        double Vsum = 0;
        double maxval = -1;
        unsigned int maxbin;
        for(unsigned int j = 0; j < n; j++) {
            Vsum += V[vector][j];
            if(V[vector][j] > maxval) {
                maxbin = j;
                maxval = V[vector][j];
            }
        }
        double cumulative_correction = 0;
        for(unsigned int j = 0; j < maxbin; j++) {
            const double correction = (V[vector][j] - V[vector][j+1]);
            if(correction > 0) {
                V[vector][j+1] += correction;
                cumulative_correction += correction;
            }
        }
        for(unsigned int j = n-1; j > maxbin; j--) {
            const double correction = (V[vector][j] - V[vector][j-1]);
            if(correction > 0) {
                V[vector][j-1] += correction;
                cumulative_correction += correction;
            }
        }
        const double scale = Vsum/(Vsum+cumulative_correction);
        for(unsigned int j = 0; j < n; j++) {
            V[vector][j] *= scale;
        }

    }
    P->SetV(V);
}

void Updating::ForceUnimodalY(Pidrix *P) {
    const unsigned int m = P->Rows();
    const unsigned int rank = P->Rank();

    TMatrixD U = P->GetU();

    for(int vector = 0; vector < rank; vector++) {

        double Usum = 0;
        double maxval = -1;
        unsigned int maxbin;
        for(unsigned int i = 0; i < m; i++) {
            Usum += U[i][vector];
            if(U[i][vector] > maxval) {
                maxbin = i;
                maxval = U[i][vector];
            }
        }
        double cumulative_correction = 0;
        for(unsigned int i = 0; i < maxbin; i++) {
            const double correction = (U[i][vector] - U[i+1][vector]);
            if(correction > 0) {
                U[i+1][vector] += correction;
                cumulative_correction += correction;
            }
        }
        for(unsigned int i = m-1; i > maxbin; i--) {
            const double correction = (U[i][vector] - U[i-1][vector]);
            if(correction > 0) {
                U[i-1][vector] += correction;
                cumulative_correction += correction;
            }
        }

        const double scale = Usum/(Usum+cumulative_correction);
        for(unsigned int i = 0; i < m; i++) {
            U[i][vector] *= scale;
        }
    }
    P->SetU(U);
}

void Updating::ForceUnimodal(Pidrix *P) {
    ForceUnimodalX(P);
    ForceUnimodalY(P);
}

void Updating::SetYield(Pidrix *P, unsigned int vector, double yield) {
    TMatrixD U = P->GetU();
    const unsigned int m = P->Rows();
    TMatrixD Ucolumn(m,1);

    TMatrixDColumn(Ucolumn, 0) = TMatrixDColumn(U, vector);
    const double old_yield = Quantifying::Yield(P, vector);
    const double scale = yield/old_yield;

    for(unsigned int i = 0; i < m; i++) {
        U[i][vector] *= scale;
    }
    P->SetU(U);
    //Normalize(P);
}

void Updating::RandomizeAmplitudes(Pidrix *P, double randomness) {
    if(randomness > 1) { randomness = 1; }
    if(randomness < 0) { randomness = 0; }

    const unsigned int rank = P->Rank();
    const double integral = P->Integral();
    TMatrixD U = P->GetU();
    TMatrixD V = P->GetV();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    TMatrixD Ucolumn(m,1), Vrow(1,n);


    for(int vector = 0; vector < rank; vector++) {
        TMatrixDColumn(Ucolumn, 0) = TMatrixDColumn(U, vector);
        TMatrixDRow(Vrow, 0) = TMatrixDRow(V, vector);
        const double yield = (Ucolumn*Vrow).Sum();
        const double new_yield = (gRandom->Uniform(0,yield*2.0))*randomness 
                                + (1.0-randomness)*yield;
        const double scale = new_yield/yield;
        for(unsigned int i = 0; i < m; i++) {
            U[i][vector] *= scale;
        }
    }
    P->SetU(U);
    P->SetV(V);
    Normalize(P);
}

double Updating::GibbsSample(Pidrix *P, unsigned int iterations, double uncertainty_scale, bool rootN_scaling) {
    TMatrixD U = P->GetU();
    TMatrixD V = P->GetV();
    const TMatrixD T = P->GetT();
    const unsigned int m = P->Rows();
    const unsigned int n = P->Columns();
    const unsigned int rank = P->Rank();

    double *Usums = new double [rank];
    double *Vsums = new double [rank];
    for(unsigned int vector = 0; vector < rank; vector++) {
        Usums[vector] = 0;
        Vsums[vector] = 0;
        for(unsigned int i = 0; i < m; i++) {
            Usums[vector] += U[i][vector];
        }
        for(unsigned int j = 0; j < m; j++) {
            Vsums[vector] += V[vector][j];
        }
    }

    //build a log likelihood matrix for each bin
    TMatrixD A = U*V;
    //the m row and n column will hold the sums
    TMatrixD LL(m+1, n+1);
    const double log_zero = log(0);
    for(unsigned int i = 0; i < m; i++) {
        double row_sum = 0;
        for(unsigned int j = 0; j < n; j++) {
            LL[i][j] = log(TMath::Poisson(T[i][j], A[i][j]));
            if(LL[i][j] != log_zero) {
                row_sum += LL[i][j];
            }
        }
        LL[i][n] = row_sum;
    }
    for(unsigned int j = 0; j < n; j++) {
        double column_sum = 0;
        for(unsigned int i = 0; i < m; i++) {
            if(LL[i][j] != log_zero) {
                column_sum += LL[i][j];
            }
        }
        LL[m][j] = column_sum;
    }

    unsigned int successes = 0, failures = 0;

    for(unsigned int iteration = 0; iteration < iterations; iteration++) {
        for(unsigned int vector = 0; vector < rank; vector++) {
            for(unsigned int i = 0; i < m; i++) {
                //interpolation is done so that the transition kernel is symmetric
                const double left_neighbor = i>0 ? U[i-1][vector] : 0;
                const double right_neighbor = i+1<m ? U[i+1][vector] : 0;
                const double interpolated = 0.5*(left_neighbor + right_neighbor);
                // +1 so it can't get stuck at 0
                const double fractional_sigma = uncertainty_scale*0.1*(rootN_scaling?1.0/sqrt(interpolated*Vsums[vector] + 1.0):1.0);
                const double sigma = fractional_sigma*(interpolated + 1.0/Vsums[vector]);
                const double delta = P->RandomGaussian(0, sigma);
                const double new_value = U[i][vector] + delta;
                //negative steps are invalid
                if(new_value < 0) {
                    continue;
                }

                double log_likelihood = 0;
                double *log_likelihood_cache = new double [n];
                for(unsigned int j = 0; j < n; j++) {
                    const double new_approximation = A[i][j] + delta*V[vector][j];
                    log_likelihood_cache[j] = log(TMath::Poisson(T[i][j], new_approximation));
                    if(log_likelihood_cache[j] != log_zero) {
                        log_likelihood += log_likelihood_cache[j];
                    }
                    else if(LL[i][j] != log_zero) {
                        log_likelihood += LL[i][j];
                    }
                }

                const double alpha = exp(log_likelihood - LL[i][n]);
                if(alpha > 1 || alpha > P->RandomUniform(0, 1)) {
                    U[i][vector] = new_value;
                    Usums[vector] += delta;
                    for(unsigned int j = 0; j < n; j++) {
                        A[i][j] += delta*V[vector][j];
                        LL[i][j] = log_likelihood_cache[j];
                    }
                    LL[i][n] = log_likelihood;
                    successes++;
                }
                else {
                    failures++;
                }
                delete [] log_likelihood_cache;
            }
            for(unsigned int j = 0; j < n; j++) {
                //interpolation is done so that the transition kernel is symmetric
                const double left_neighbor = j>0 ? V[vector][j-1] : 0;
                const double right_neighbor = j+1<n ? V[vector][j+1] : 0;
                const double interpolated = 0.5*(left_neighbor + right_neighbor);
                // +1 so it can't get stuck at 0
                const double fractional_sigma = uncertainty_scale*0.1*(rootN_scaling?1.0/sqrt(interpolated*Usums[vector] + 1.0):1.0);
                const double sigma = fractional_sigma*(interpolated + 1.0/Usums[vector]);
                const double delta = P->RandomGaussian(0, sigma);
                const double new_value = V[vector][j] + delta;
                //negative steps are invalid
                if(new_value < 0) {
                    continue;
                }

                double log_likelihood = 0;
                double *log_likelihood_cache = new double [m];
                for(unsigned int i = 0; i < m; i++) {
                    const double new_approximation = A[i][j] + delta*U[i][vector];
                    log_likelihood_cache[i] = log(TMath::Poisson(T[i][j], new_approximation));
                    if(log_likelihood_cache[i] != log_zero) {
                        log_likelihood += log_likelihood_cache[i];
                    }
                    else if(LL[i][j] != log_zero) {
                        log_likelihood += LL[i][j];
                    }
                }

                const double alpha = exp(log_likelihood - LL[m][j]);
                if(alpha > 1 || alpha > P->RandomUniform(0, 1)) {
                    V[vector][j] = new_value;
                    Vsums[vector] += delta;
                    for(unsigned int i = 0; i < m; i++) {
                        A[i][j] += delta*U[i][vector];
                        LL[i][j] = log_likelihood_cache[i];
                    }
                    LL[m][j] = log_likelihood;
                    successes++;
                }
                else {
                    failures++;
                }
                delete [] log_likelihood_cache;
            }
        }
    }

    P->SetU(U);
    P->SetV(V);
    delete [] Usums;
    delete [] Vsums;

    return double(successes)/double(successes+failures);
}
