#include "ExampleFunction.h"
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;
// minimize 3*x^2 + 2*x*y + 2*y^2 + 7

ExampleFunction::ExampleFunction(Placement &placement)
    :_placement(placement)
{
    // for WA
    _compute = false;
    _gamma = (_placement.boundryTop() - _placement.boundryBottom()) / 600;
    _wl = 0;
    _density = 0;
    // positive_weight_x.resize(_placement.numNets(),0);
    // positive_x.resize(_placement.numNets(),0);
    // negative_weight_x.resize(_placement.numNets(),0);
    // negative_x.resize(_placement.numNets(),0);
    // positive_weight_y.resize(_placement.numNets(),0);
    // positive_y.resize(_placement.numNets(),0);
    // negative_weight_y.resize(_placement.numNets(),0);
    // negative_y.resize(_placement.numNets(),0);

    //for density
    _computeD = false;
    _edgenum = 14;
    _binW = (_placement.boundryRight() - _placement.boundryLeft()) / _edgenum;
    _binH = (_placement.boundryTop() - _placement.boundryBottom()) / _edgenum;
    _binM = 0;
    for(unsigned i = 0;i < _placement.numModules();i++) _binM += _placement.module(i).area();
    _binM /= ((_placement.boundryRight() - _placement.boundryLeft()) * (_placement.boundryTop() - _placement.boundryBottom()));
    binD.resize(pow(_edgenum,2),0);

    _lamda = 0;
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
    // f = 3*x[0]*x[0] + 2*x[0]*x[1] + 2*x[1]*x[1] + 7; // objective function
    // g[0] = 6*x[0] + 2*x[1]; // gradient function of X
    // g[1] = 2*x[0] + 4*x[1]; // gradient function of Y
    g.clear();
    g.resize(2*_placement.numModules(),0);
    _density = 0;
    // f = compute_WA_g(x,g) + compute_D_g(x,g);
    // compute_WA_g(x,g); // gradient of wirelength model
    // compute_D_g(x,g);

    vector<double> temp;
    temp.resize(4*_placement.numModules(),0);
    f = 0.0;
    // LSE
    for(unsigned i = 0;i < _placement.numModules();i++)
    {
        // cout << "x" << x[2*i] << ", " << x[2*i+1] << endl;
        temp[4*i] = exp(x[2*i]/_gamma);
        temp[4*i + 1] = exp((-1)*x[2*i]/_gamma);
        temp[4*i + 2] = exp(x[2*i+1]/_gamma);
        temp[4*i + 3] = exp((-1)*x[2*i+1]/_gamma);
        // cout << "temp" << temp[4*i] << ", " << temp[4*i+1]<<endl;
    }
    for(unsigned i = 0;i < _placement.numNets();i++)
    {
        double L1 = 0;
        double L2 = 0;
        double L3 = 0;
        double L4 = 0;
        if(_placement.net(i).numPins() == 0) continue;
        for(unsigned j = 0;j < _placement.net(i).numPins();j++)
        {
            unsigned id = _placement.net(i).pin(j).moduleId();
            L1 += temp[4*id];
            L2 += temp[4*id + 1];
            L3 += temp[4*id + 2];
            L4 += temp[4*id + 3];
        }
        // cout << L1 << " " << L2 << endl;
        f += log(L1) + log(L2) + log(L3) + log(L4);

        for(unsigned j = 0;j < _placement.net(i).numPins();j++)
        {
            unsigned id = _placement.net(i).pin(j).moduleId();
            g[2*id] += temp[4*id] / (_gamma*L1);  
            g[2*id] -= temp[4*id+1] / (_gamma*L2);  
            g[2*id+1] += temp[4*id+2] / (_gamma*L3);  
            g[2*id+1] -= temp[4*id+3] / (_gamma*L4);    
        }
    }
    // f *= _gamma;
    // cout << "WL2 = "<<f<<endl;
    if(_lamda == 0) return;

    // cout << "start compute_D_g\n";
    double _ax,_bx,_ay,_by;
    double _dx,_dy,_ratio;
    double _px,_py;
    double _orix,_oriy;
    double _cx,_cy;
    double _bincx,_bincy;
    double counter = 0;
    binD.clear();
    binD.resize(pow(_edgenum,2),0);
    vector<double> tempg;
    tempg.resize(g.size(),0);

    for(unsigned i = 0;i < _edgenum;i++)
    {
        for(unsigned j = 0;j < _edgenum;j++)
        {
            for(unsigned k = 0;k < _placement.numModules();k++)
            {
                Module m = _placement.module(k);

                _ax = 4 / ( (m.width()+2*_binW) * (m.width()+4*_binW) );
                _bx = 2 / (_binW * (m.width()+4*_binW));
                _ay = 4 / ( (m.height()+2*_binH) * (m.height()+4*_binH) );
                _by = 2 / (_binH * (m.height()+4*_binH));

                _ratio = m.area() / (_binW*_binH);

                _cx = x[2*k] + m.width()/2;
                _cy = x[2*k + 1] + m.height()/2;
                _bincx = (i+0.5)*_binW + _placement.boundryLeft();
                _bincy = (j+0.5)*_binH + _placement.boundryBottom();

                if(_bincx <= _cx) _dx = _cx - (_bincx);
                else _dx = _bincx - _cx;
                if(_bincy <= _cy) _dy = _cy - (_bincy);
                else _dy = _bincy - _cy;


                // _dx = x[2*k] - ((i+0.5)*_binW - _placement.boundryLeft());
                // _dy = x[2*k+1] - ((j+0.5)*_binH - _placement.boundryBottom());
                // _orix = _dx;
                // _oriy = _dy;
                // _dx = abs(_dx);
                // _dy = abs(_dy);
                // if(_dx < 0 || _dy < 0) cout << "minus\n";

                if(_dx <= (m.width()/2 + _binW))    
                    _px = 1 - (_ax * pow(_dx,2));
                else if (_dx <= (m.width()/2 + 2*_binW))
                    _px = _bx * pow((_dx - m.width()/2 - 2*_binW),2);
                else
                    _px = 0;

                if(_dy <= (m.height()/2 + _binH))    
                    _py = 1 - (_ay * pow(_dy,2));
                else if (_dy <= (m.height()/2 + 2*_binH))
                    _py = _by * pow((_dy - m.height()/2 - 2*_binH),2);
                else
                    _py = 0;

                if(_dx <= (m.width()/2 + _binW))
                {
                    // cout << "1\n";
                    tempg[2*k] = _ratio * (-2) * _ax * _dx * _py;
                    if(_bincx > _cx) tempg[2*k] *= -1;
                }
                else if (_dx <= (m.width()/2 + 2*_binW))
                {
                    // cout << "2\n";
                    tempg[2*k] = _ratio * 2 * _bx * (_dx - m.width()/2 - 2*_binW) * _py;
                    if(_bincx > _cx) tempg[2*k] *= -1;
                }
                else
                    tempg[2*k] = 0;
                
                if(_dy <= (m.height()/2 + _binH))
                {
                    // cout << "3\n";
                    tempg[2*k+1] = _ratio * (-2) * _ay * _dy * _px;
                    if(_bincy > _cy) tempg[2*k+1] *= -1;
                }
                else if (_dy <= (m.height()/2 + 2*_binH))
                {
                    // cout << "4\n";
                    tempg[2*k+1] = _ratio * 2 * _by * (_dy - m.height()/2 - 2*_binH) * _px;
                    if(_bincy > _cy) tempg[2*k+1] *= -1;
                }
                else
                    tempg[2*k+1] = 0;
                // if(_px != 0 || _py != 0)
                //     cout << "px = "<<_px<<", py = "<<_py<<endl;
                binD[i+_edgenum*j] += _ratio*_px*_py;
            }
            _density += _lamda * pow((binD[i+_edgenum*j] - _binM),2);
            if(binD[i+_edgenum*j] <= _binM) counter++;
            for(unsigned a = 0;a < _placement.numModules();a++)
            {
                g[2*a] += _lamda * 2 * (binD[i+_edgenum*j] - _binM) * tempg[2*a]; 
                g[2*a + 1] += _lamda * 2 * (binD[i+_edgenum*j] - _binM) * tempg[2*a + 1]; 
            }
        }
    }
    // _density *= _lamda;
    // cout << "density = "<<_density<<endl;
    // for(int i = 0;i < _edgenum;i++)
    // {
    //     for(int j = 0;j < _edgenum;j++)
    //     {
    //         cout << setprecision(3) << binD[i+_edgenum*j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "counter = " << counter << ", ratio = " << counter / pow(_edgenum,2) << endl; 
    _wl = counter;
    f += _density;
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{
    // f = 3*x[0]*x[0] + 2*x[0]*x[1] + 2*x[1]*x[1] + 7; // objective function
    // f = compute_WA(x) + compute_D(x); // wirelength model

    vector<double> temp;
    temp.resize(4*_placement.numModules(),0);
    f = 0.0;
    // LSE
    for(unsigned i = 0;i < _placement.numModules();i++)
    {
        temp[4*i] = exp(x[2*i]/_gamma);
        temp[4*i + 1] = exp((-1)*x[2*i]/_gamma);
        temp[4*i + 2] = exp(x[2*i+1]/_gamma);
        temp[4*i + 3] = exp((-1)*x[2*i+1]/_gamma);
    }
    for(unsigned i = 0;i < _placement.numNets();i++)
    {
        double L1 = 0;
        double L2 = 0;
        double L3 = 0;
        double L4 = 0;
        if(_placement.net(i).numPins() == 0) continue;
        for(unsigned j = 0;j < _placement.net(i).numPins();j++)
        {
            unsigned id = _placement.net(i).pin(j).moduleId();
            L1 += temp[4*id];
            L2 += temp[4*id + 1];
            L3 += temp[4*id + 2];
            L4 += temp[4*id + 3];
        }
        f += log(L1) + log(L2) + log(L3) + log(L4);
    }
    // f *= _gamma;
    
    if(_lamda == 0) return;

    double _ax,_bx,_ay,_by;
    double _dx,_dy,_ratio;
    double _px,_py;
    double _cx,_cy;
    double _bincx,_bincy;
    double counter = 0;
    _density = 0;
    binD.clear();
    binD.resize(pow(_edgenum,2),0);

    for(unsigned i = 0;i < _edgenum;i++)
    {
        for(unsigned j = 0;j < _edgenum;j++)
        {
            for(unsigned k = 0;k < _placement.numModules();k++)
            {
                Module m = _placement.module(k);

                _ax = 4 / ( (m.width()+2*_binW) * (m.width()+4*_binW) );
                _bx = 2 / (_binW * (m.width()+4*_binW));
                _ay = 4 / ( (m.height()+2*_binH) * (m.height()+4*_binH) );
                _by = 2 / (_binH * (m.height()+4*_binH));

                _ratio = m.area() / (_binW*_binH);

                _cx = x[2*k] + m.width()/2;
                _cy = x[2*k + 1] + m.height()/2;
                _bincx = (i+0.5)*_binW + _placement.boundryLeft();
                _bincy = (j+0.5)*_binH + _placement.boundryBottom();

                if(_bincx <= _cx) _dx = _cx - (_bincx);
                else _dx = _bincx - _cx;
                if(_bincy <= _cy) _dy = _cy - (_bincy);
                else _dy = _bincy - _cy;

                // _dx = x[2*k] - ((i+0.5)*_binW - _placement.boundryLeft());
                // _dy = x[2*k+1] - ((j+0.5)*_binH - _placement.boundryBottom());
                // _orix = _dx;
                // _oriy = _dy;
                // _dx = abs(_dx);
                // _dy = abs(_dy);

                if(_dx <= (m.width()/2 + _binW))    
                    _px = 1 - (_ax * pow(_dx,2));
                else if (_dx <= (m.width()/2 + 2*_binW))
                    _px = _bx * pow((_dx - m.width()/2 - 2*_binW),2);
                else
                    _px = 0;

                if(_dy <= (m.height()/2 + _binH))    
                    _py = 1 - (_ay * pow(_dy,2));
                else if (_dy <= (m.height()/2 + 2*_binH))
                    _py = _by * pow((_dy - m.height()/2 - 2*_binH),2);
                else
                    _py = 0;

                binD[i+_edgenum*j] += _ratio*_px*_py;
            }
            _density += _lamda * pow((binD[i+_edgenum*j] - _binM),2);
            if(binD[i+_edgenum*j] <= _binM) counter++;
        }
    }
    // _computeD = true;
    // cout << "density = "<<_density<<endl;
    // _density *= _lamda;
// cout << "density = "<<_density<<endl;
    // for(int i = 0;i < _edgenum;i++)
    // {
    //     for(int j = 0;j < _edgenum;j++)
    //     {
    //         cout << setprecision(3) << binD[i+_edgenum*j] << " ";
    //     }
    //     cout << endl;
    // }
    f += _density;
    // cout << "counter = " << counter << ", ratio = " << counter / pow(_edgenum,2) << endl; 
    _wl = counter;
}

unsigned ExampleFunction::dimension()
{
    return _placement.numModules()*2; // num_blocks*2 
    // each two dimension represent the X and Y dimensions of each block
}

void ExampleFunction::debug()
{
    for(int i = 0;i < _edgenum;i++)
    {
        for(int j = 0;j < _edgenum;j++)
        {
            cout << setprecision(3) << binD[i+_edgenum*j] << " ";
        }
        cout << endl;
    }
    cout << "counter = " << _wl << ", ratio = " << _wl/pow(_edgenum,2) << endl;
}