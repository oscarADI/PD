#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Placement.h"

class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(Placement &placement);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();

    double compute_WA(const vector<double> &x);
    double compute_WA_g(const vector<double> &x,vector<double> &g);
    double compute_D(const vector<double> &x);
    double compute_D_g(const vector<double> &x,vector<double> &g);
    void setlamda(double i) {_lamda = i;}

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
