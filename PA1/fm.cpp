#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

void parser (const char* argv, double& degree, vector< vector<string> >& Nets)
{
    ifstream infile;
    infile.open(argv,ifstream::in);
    if(!infile.is_open())
    {
        cerr << "Input file doesn't exist !" << endl;
        return;
    }
    //parser(infile);
    string line;
    string net;
    vector<string> temp;
    while ( getline(infile,line))
    {
        // cout << line << endl;
        if(line[0] == '0')
        {
            degree = stod(line);
            continue;
        }
        size_t first = 0;
        size_t second = 0;
        if(line[0] == 'N')
        {
            // for(vector<string>::iterator iter = temp.begin(); iter != temp.end(); iter++)
            // {
            //     cout << *iter << endl;
            // }
            temp.clear();
        }
        while(true)
        {
            first = line.find('c', first+1);
            if(first == string::npos) break;
            second = line.find(' ',first+1);
            net = line.substr(first,second-first);
            temp.push_back(net);
        }       
    }

}

int main(int argc, char* argv[])
{
    double degree;
    vector< vector<string> > Nets;
    parser(argv[1],degree,Nets);

}
