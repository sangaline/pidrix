#include "initialization.h"
#include "plotting.h"

#include "Pidrix.h"

#include "TGraph.h"
#include "TF1.h"

#include "TVectorD.h"
#include "TMatrixD.h"

#include "math.h"

using namespace Initialization;

//Vector initializations
//======================

double Initialization::UniformRandomVectors(Pidrix *P) {
    //if each vector sums to this then our integral should match
    const double vector_sum = sqrt(P->Integral()/double(P->Rank()));
    const double Uvector_max_entry = 2.0*vector_sum/double(P->Rows());
    const double Vvector_max_entry = 2.0*vector_sum/double(P->Columns());

    double integral = 0;
    for(unsigned int vector = 0; vector < P->Rank(); vector++) {

        double Usum = 0, Vsum = 0;
        for(unsigned int i = 0; i < P->Rows(); i++) {
            P->SetU(i, vector, P->RandomUniform(0, Uvector_max_entry));
            Usum += P->GetU(i, vector);
        }
        for(unsigned int j = 0; j < P->Columns(); j++) {
            P->SetV(vector, j, P->RandomUniform(0, Vvector_max_entry));
            Vsum += P->GetV(vector, j);
        }

        integral += Usum*Vsum;
    }

    return integral;
}


//Rank initializations
//====================

unsigned int Initialization::SVDThresholdRank(Pidrix *P, double threshold) {
    //Guess the rank
    const TVectorD S = P->SVDSigma();
    const int Slen = S.GetNoElements();

    double absolute_total = 0;
    for(unsigned int i = 0; i < Slen; i++) {
        absolute_total += S[i];
    }
    unsigned int rank = 1;
    double current_total = 0;
    for(unsigned int i = 0; i < Slen; i++) {
        current_total += S[i];
        if(current_total / absolute_total > threshold) break;
        rank++;
    }
    return P->SetRank(rank);
}

unsigned int Initialization::SVDLinearIntersectionRank(Pidrix *P) {
    TGraph *sv = Plotting::SVGraph(P);
    unsigned int Slen = sv->GetN();

    TF1 *line = new TF1("line", "[0]+[1]*x", 0, Slen + 1.0);
    sv->Fit(line, "QN");

    unsigned int rank = 0;
    for(int i = 0; i < Slen; i++) {
        double x, y;
        sv->GetPoint(i, x, y);
        if( y < line->Eval(i+1) ) break;
        rank++;
    }

    return P->SetRank(rank);
}
