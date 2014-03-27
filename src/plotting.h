#ifndef plotting_h
#define plotting_h

class Pidrix;
class TGraph;
class TH2D;
class TH1D;

namespace Plotting {
    TGraph* SVGraph(const Pidrix* P, TGraph *t = 0);
    TH2D* Approximation(const Pidrix* P, TH2D* h = 0, const char* name = "approximation");
    TH2D* Target(const Pidrix* P, TH2D* h = 0, const char* name = "target");
    TH1D* DistributionX(const Pidrix* P, unsigned int vector, TH1D* h = 0);
    TH1D* DistributionY(const Pidrix* P, unsigned int vector, TH1D* h = 0);
    TH2D* DistributionXY(const Pidrix* P, unsigned int vector, TH2D* h = 0);

};

#endif
