#ifndef plotting_h
#define plotting_h

class Pidrix;
class TGraph;
class TH2D;

namespace Plotting {
    TGraph* SVGraph(Pidrix* P, TGraph *t = 0);
    TH2D* Approximation(Pidrix* P, TH2D* h = 0, const char* name = "approximation");

};

#endif
