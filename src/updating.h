#ifndef updating_h
#define updating_h

class Pidrix;

namespace Updating {
    void MultiplicativeEuclidian(Pidrix *P, const unsigned int iterations = 1);
    void MultiplicativeKL(Pidrix *P, const unsigned int iterations = 1);

    void Normalize(Pidrix *P);
    void ScaleX(Pidrix *P, double factor);
    void ScaleY(Pidrix *P, double factor);
    void Scale(Pidrix *P, double factor);
    void AddNoiseX(Pidrix *P, double fraction = 0.01);
    void AddNoiseY(Pidrix *P, double fraction = 0.01);
    void AddNoise(Pidrix *P, double fraction = 0.01);
    void Smear(Pidrix *P, unsigned int iterations = 1, double neighbor_fraction = 0.25);
    void SmearX(Pidrix *P, unsigned int iterations = 1, double neighbor_fraction = 0.25);
    void SmearY(Pidrix *P, unsigned int iterations = 1, double neighbor_fraction = 0.25);

    void ForceGaussianX(Pidrix *P);
    void ForceGaussianY(Pidrix *P);
    void ForceGaussian(Pidrix *P);
};

#endif
