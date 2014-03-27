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
};

#endif
