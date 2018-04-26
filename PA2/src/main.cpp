#include<iostream>
#include"floorplan.h"

using namespace std;

int main(int argc, char* argv[])
{
    Floorplan test;
    test.seta(argv[1]);
    test.parser_block(argv[2]);
    test.parser_net(argv[3]);
    // test.debug();
    test.init_tree();
    test.packing();
    test.compute_norm();
    test.SA();
}
