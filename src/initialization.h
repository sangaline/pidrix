#ifndef initialization_h
#define initialization_h
class Pidrix;

namespace Initialization {
    double UniformRandomVectors(Pidrix *P, bool randomize_scale = true);

    unsigned int SVDThresholdRank(Pidrix *P, double threshold = 0.995);
    unsigned int SVDLinearIntersectionRank(Pidrix *P);
};
#endif
