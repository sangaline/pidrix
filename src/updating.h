#ifndef updating_h
#define updating_h

class Pidrix;

namespace Updating {
    void MultiplicativeEuclidian(Pidrix *P, const unsigned int iterations = 1);
    void MultiplicativeKL(Pidrix *P, const unsigned int iterations = 1);

    void Normalize(Pidrix *P);
};

#endif
