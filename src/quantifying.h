#ifndef quantifying_h
#define quantifying_h

class Pidrix;

namespace Quantifying {
    double Yield(Pidrix *P, unsigned int vector);
    double MeanX(Pidrix *P, unsigned int vector);
    double MeanY(Pidrix *P, unsigned int vector);
    double VarianceX(Pidrix *P, unsigned int vector);
    double VarianceY(Pidrix *P, unsigned int vector);
    double StandardDeviationX(Pidrix *P, unsigned int vector);
    double StandardDeviationY(Pidrix *P, unsigned int vector);
    double X(Pidrix *P, unsigned int j);
    double Y(Pidrix *P, unsigned int i);
};

#endif
