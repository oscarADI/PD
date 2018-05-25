#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Placement.h"
#define _GLIBCXX_USE_CXX11_ABI 0

class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(Placement &placement);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();

    void setlamda(double i) {_lamda = i;}
    void debug();

private:
    Placement& _placement;
    bool       _compute;
    double     _gamma;
    double     _wl;
    vector<double> positive_weight_x;
    vector<double> positive_x;
    vector<double> negative_weight_x;
    vector<double> negative_x;
    vector<double> positive_weight_y;
    vector<double> positive_y;
    vector<double> negative_weight_y;
    vector<double> negative_y;

    double     _edgenum;
    double     _binW,_binH;
    double     _binM;
    double     _density;
    bool       _computeD;
    vector<double> binD;

    double     _lamda;
};
#endif // EXAMPLEFUNCTION_H
