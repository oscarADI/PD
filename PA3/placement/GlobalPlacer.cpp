#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"

GlobalPlacer::GlobalPlacer(Placement &placement)
	:_placement(placement)
{

}


void GlobalPlacer::place()
{
	///////////////////////////////////////////////////////////////////
	// The following example is only for analytical methods.
	// if you use other methods, you can skip and delete it directly.
	//////////////////////////////////////////////////////////////////

	ExampleFunction ef(_placement); // require to define the object function and gradient function

    vector<double> x; // solution vector, size: num_blocks*2 
                         // each 2 variables represent the X and Y dimensions of a block
    // x[0] = 100; // initialize the solution vector
    // x[1] = 100;
    x.resize(2*_placement.numModules(),0);
    // for(unsigned i = 0;i < _placement.numModules();i++)
    // {
    //     x[2*i] = ( _placement.boundryRight() + _placement.boundryLeft() )/ 2;
    //     x[2*i + 1] = (_placement.boundryTop() + _placement.boundryBottom()) / 2;
    //     x[2*i] = 0;
    //     x[2*i + 1] = 0;
    //     _placement.module(i).setPosition(x[2*i],x[2*i+1]);
    // }

    NumericalOptimizer no(ef);
    
    for(unsigned i = 0;i < 3;i++)
    {
        cout << "i = " << i << endl;
        ef.setlamda(i*500);
        no.setX(x); // set initial solution
        no.setNumIteration(35); // user-specified parameter
        no.setStepSizeBound(( _placement.boundryRight() - _placement.boundryLeft() )*5); // user-specified parameter
        no.solve(); // Conjugate Gradient solver

        for(unsigned j = 0;j < _placement.numModules();j++)
        {
            double newx = no.x(2*j);
            double newy = no.x(2*j+1);

            if(newx < _placement.boundryLeft()) 
                newx = _placement.boundryLeft();
            else if(newx + _placement.module(j).width() > _placement.boundryRight())
                newx = _placement.boundryRight() - _placement.module(j).width();
            
            if(newy < _placement.boundryBottom()) 
                newy = _placement.boundryBottom();
            else if(newy + _placement.module(j).height() > _placement.boundryTop())
                newy = _placement.boundryTop() - _placement.module(j).height();

            _placement.module(j).setPosition(newx,newy);
            x[2*j] = newx;
            x[2*j + 1] = newy;
        }
        if(i == 0)  plotPlacementResult( "init0.plt" );
        else if(i == 1) plotPlacementResult( "init1.plt" );
        else if(i == 2) plotPlacementResult( "init2.plt" );
        else if(i == 3) plotPlacementResult( "init3.plt" );
        else if(i == 4) plotPlacementResult( "init4.plt" );
    }

    cout << "Current solution:" << endl;
    for (unsigned i = 0; i < no.dimension(); i++) {
        cout << "x[" << i << "] = " << no.x(i) << endl;
    }
    cout << "Objective: " << no.objective() << endl;
	////////////////////////////////////////////////////////////////


}


void GlobalPlacer::plotPlacementResult( const string outfilename, bool isPrompt )
{
    ofstream outfile( outfilename.c_str() , ios::out );
    outfile << " " << endl;
    outfile << "set title \"wirelength = " << _placement.computeHpwl() << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl << endl;
    outfile << "# bounding box" << endl;
    plotBoxPLT( outfile, _placement.boundryLeft(), _placement.boundryBottom(), _placement.boundryRight(), _placement.boundryTop() );
    outfile << "EOF" << endl;
    outfile << "# modules" << endl << "0.00, 0.00" << endl << endl;
    for( size_t i = 0; i < _placement.numModules(); ++i ){
        Module &module = _placement.module(i);
        plotBoxPLT( outfile, module.x(), module.y(), module.x() + module.width(), module.y() + module.height() );
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if( isPrompt ){
        char cmd[ 200 ];
        sprintf( cmd, "gnuplot %s", outfilename.c_str() );
        if( !system( cmd ) ) { cout << "Fail to execute: \"" << cmd << "\"." << endl; }
    }
}

void GlobalPlacer::plotBoxPLT( ofstream& stream, double x1, double y1, double x2, double y2 )
{
    stream << x1 << ", " << y1 << endl << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl << endl;
}
