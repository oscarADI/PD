#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <assert.h>
using namespace std;

class Net {
    public:
        Net() 
        {
            index = 0;
        }
        Net(int n)
        {
            index = n;
        }
        ~Net() {}

        // initialize
        void pushcell(int n) { cells.push_back(n); }
        int cellsize() { return cells.size(); }
        int cell(int n) { return cells[n]; }
        int netindex() { return index; }
        int partsize(int i) { if(!i) return cellsin0.size(); else return cellsin1.size(); }
        int partcellbegin(int i) { if(!i) return cellsin0[0]; else return cellsin1[0]; }
        void pushpart(int n, int i) 
        {
            if(!i) cellsin0.push_back(n);
            else   cellsin1.push_back(n);
        }
        void movepart(int n, int i)
        {
            if(!i)
            {
                vector<int>::iterator it;
                it = find(cellsin0.begin(),cellsin0.end(),n);
                cellsin0.erase(it);
                cellsin1.push_back(n);
            }
            else
            {
                vector<int>::iterator it;
                it = find(cellsin1.begin(),cellsin1.end(),n);
                cellsin1.erase(it);
                cellsin0.push_back(n);
            }
        }

        // debug
        void print() {
            cout << "n" << index << " "; 
            for(vector<int>::iterator iter = cells.begin(); iter != cells.end(); iter++)
                cout << "c" << *iter << " ";
            cout<<endl;
        }

    private:
        int             index;
        vector<int>     cells;
        vector<int>     cellsin0;
        vector<int>     cellsin1;
};

class Cell {
    public:
        Cell()
        {
            index = locked = gain = 0;
            partition = -1;
        }
        Cell(int n)
        {
            index = n;
            locked = gain = 0;
            partition = partition_init = -1;
        }
        ~Cell() {}

        // initialize
        void pushnet(int n) { nets.push_back(n); }
        int Psize() { return nets.size(); }
        int cellindex() { return index; }
        void assign_part(int n) { partition = partition_init = n; }
        int gain_num() { return gain; }
        void change_gain(int n) { gain += n; }
        int part_num() { return partition; }
        int net(int n) { return nets[n]; }
        int lock() {return locked;}
        void makelock() {locked = 1;}
        void change_part() {
            if(partition) partition = 0; 
            else partition = 1;
            locked = 1;
        }

        // debug
        void print() {
            cout << "c" << index << " "; 
            for(vector<int>::iterator iter = nets.begin(); iter != nets.end(); iter++)
                cout << "n" << *iter << " ";
            cout<<endl;
        }

    private:
        int             index;
        vector<int>     nets;
        int             partition;
        int             partition_init;
        int             locked;
        int             gain;
};

class FM {
    public:
        FM() {
            nCells = 0;
            Pmax = 0;
            degree = lowerbound = upperbound = 0.0;
        }
        ~FM() {}

        void            parser(const char*);
        int             netsize() { return netlist.size(); }
        void            compute_gain();
        void            update_gain(Cell*);
        void            iteration();
        bool            balanced(int);

        // debug
        void            print_cells();

    private:
        int                                 nCells;
        vector<Net*>                        netlist;
        vector<Cell*>                       celllist;
        double                              degree;
        vector< vector<int> >               list0;
        vector< vector<int> >               list1;
        int                                 Pmax;
        double                              lowerbound;
        double                              upperbound;
        int                                 Gmax0;
        int                                 Gmax1;
};
bool FM::balanced(int n)
{
    int size0 = list0.size();
    int size1 = list1.size();
    if(n) 
    {
        size0+=1;
        size1-=1;
    }
    else
    {
        size0-=1;
        size1+=1;
    }
    return ((size0>=lowerbound)&&(size0<=upperbound)&&(size1>=lowerbound)&&(size1<=upperbound));
}
void FM::print_cells()
{
    for(vector<Cell*>::iterator iter = celllist.begin(); iter != celllist.end(); iter++)
    {
        if((*iter))
        {
            cout << (*iter)->cellindex() << "  gain = " << (*iter)->gain_num() <<endl;
        }
    }
}

void FM::parser (const char* argv)
{
    ifstream infile;
    infile.open(argv,ifstream::in);
    if(!infile.is_open())
    {
        cerr << "Input file doesn't exist !" << endl;
        return;
    }

    string line;
    getline(infile,line);
    degree = atof(line.c_str());

    int current_net = 0;
    int cell = 0;
    while ( infile >> line )
    {
        if(line == "NET")
        {
            infile >> line;
            line = line.substr(1);
            current_net = atoi(line.c_str());
            Net* temp = new Net(current_net);

            while( infile >> line )
            {
                if(line == ";")
                {
                    netlist.push_back(temp);
                    break;
                }
                line = line.substr(1);
                cell = atoi(line.c_str());
                temp->pushcell(cell);
                nCells = max(nCells,cell);
            } 
        }   
    }
    infile.close();

    celllist.resize(nCells);
    for(vector<Net*>::iterator iter = netlist.begin(); iter != netlist.end(); iter++)
    {
        for(int i = 0 ; i < (*iter)->cellsize() ; i++)
        {
            int now = (*iter)->cell(i);
            if(celllist[now])
            {
                celllist[now]->pushnet((*iter)->netindex());
            }
            else
            {
                Cell* tmp = new Cell(now);
                tmp->pushnet((*iter)->netindex());
                celllist[now] = tmp;
            }
        }
    }
    for(vector<Cell*>::iterator iter = celllist.begin(); iter != celllist.end(); iter++)
    {
        if((*iter)) 
        {
            Pmax = max(Pmax,(*iter)->Psize());
            if((*iter)->cellindex() <= nCells/2) 
            {
                (*iter)->assign_part(0);
            }
            else 
            {
                (*iter)->assign_part(1);
            }
        }
    }

    list0.resize(2*Pmax);
    list1.resize(2*Pmax);
    Gmax0 = Gmax1 = -Pmax;
    lowerbound = (1-degree) * nCells / 2;
    upperbound = (1+degree) * nCells / 2;
}

void FM::compute_gain()
{
    for(vector<Net*>::iterator iter = netlist.begin(); iter != netlist.end(); iter++)
    {
        int p0 = 0;
        int p1 = 0;
        for(int i = 0 ; i < (*iter)->cellsize() ; i++)
        {
            Cell* tmp = celllist[(*iter)->cell(i)];
            if(tmp->part_num() == 0) 
            {
                p0++;
                (*iter)->pushpart(tmp->cellindex(),0);
            }
            else 
            {
                p1++;
                (*iter)->pushpart(tmp->cellindex(),1);
            }
        }
        if(p0 == 0)
        {
            for(int i = 0 ; i < (*iter)->cellsize() ; i++)
            {
                Cell* tmp = celllist[(*iter)->cell(i)];
                if(tmp->part_num() == 1) tmp->change_gain(-1); 
            }   
        }
        else if(p1 == 0)
        {
            for(int i = 0 ; i < (*iter)->cellsize() ; i++)
            {
                Cell* tmp = celllist[(*iter)->cell(i)];
                if(tmp->part_num() == 0) tmp->change_gain(-1); 
            }   
        }
        if(p0 == 1)
        {
            for(int i = 0 ; i < (*iter)->cellsize() ; i++)
            {
                Cell* tmp = celllist[(*iter)->cell(i)];
                if(tmp->part_num() == 1) tmp->change_gain(1); 
            }   
        }
        if(p1 == 1)
        {
            for(int i = 0 ; i < (*iter)->cellsize() ; i++)
            {
                Cell* tmp = celllist[(*iter)->cell(i)];
                if(tmp->part_num() == 0) tmp->change_gain(1); 
            }   
        }
    }
    for(vector<Cell*>::iterator iter = celllist.begin()+1; iter != celllist.end(); iter++)
    {
        if((*iter)->part_num() == 0)
        {
            list0[(*iter)->gain_num()+Pmax].push_back((*iter)->cellindex());
            Gmax0 = max(Gmax0,(*iter)->gain_num()+Pmax);
        }
        else
        {
            list1[(*iter)->gain_num()+Pmax].push_back((*iter)->cellindex());
            Gmax1 = max(Gmax1,(*iter)->gain_num()+Pmax);          
        }
    }
}

void FM::update_gain(Cell* base)
{
    int init_part = base->part_num();
    for(int i = 0 ; i < base->Psize() ; i++)
    {
        Net* tmp = netlist[base->net(i)];
        Cell* c;
        if(init_part == 0)
        {
            if(tmp->partsize(1) == 0)
                for(int j = 0; j < tmp->cellsize(); j++)
                {
                    c = celllist[tmp->cell(j)];
                    if(c->cellindex() == base->cellindex()) continue;
                    if(!c->lock())
                    {
                        int g = c->gain_num() + Pmax;
                        vector<int>::iterator it ;
                        it = find(list0[g].begin(),list0[g].end(),c->cellindex());
                        assert(it != list0[g].end());
                        list0[g].erase(it);
                        c->change_gain(1);
                        list0[g+1].push_back(c->cellindex());
                        if (g+1 > Gmax0) Gmax0 = g+1;
                    }

                }
            else if(tmp->partsize(1) == 1)
            {
                c = celllist[tmp->partcellbegin(1)];
                if(!c->lock())
                { 
                    int g = c->gain_num() + Pmax;
                    vector<int>::iterator it ;
                    it = find(list1[g].begin(),list1[g].end(),c->cellindex());
                    assert(it != list1[g].end());
                    list1[g].erase(it);
                    c->change_gain(-1);
                    list1[g-1].push_back(c->cellindex());
                    if(list1[g].size() == 0 && Gmax1 == g) Gmax1 = g-1; 
                }
            }
            tmp->movepart(base->cellindex(),0);
        }
        else
        {
            if(tmp->partsize(0) == 0)
            {
                for(int j = 0; j < tmp->cellsize(); j++)
                {
                    c = celllist[tmp->cell(j)];
                    if(c->cellindex() == base->cellindex()) continue;
                    if(!c->lock())
                    {
                        int g = c->gain_num() + Pmax;
                        vector<int>::iterator it ;
                        it = find(list1[g].begin(),list1[g].end(),c->cellindex());
                        assert(it != list1[g].end());
                        list1[g].erase(it);
                        c->change_gain(1);
                        list1[g+1].push_back(c->cellindex());
                        if (g+1 > Gmax1) Gmax1 = g+1;
                    }
                }
            }
            else if(tmp->partsize(0) == 1)
            {
                c = celllist[tmp->partcellbegin(0)];
                if(!c->lock())
                { 
                    int g = c->gain_num() + Pmax;
                    vector<int>::iterator it ;
                    it = find(list0[g].begin(),list0[g].end(),c->cellindex());
                    assert(it != list0[g].end());
                    list0[g].erase(it);
                    c->change_gain(-1);
                    list0[g-1].push_back(c->cellindex());
                    if(list0[g].size() == 0 && Gmax0 == g) Gmax0 = g-1; 
                }                
            }
            tmp->movepart(base->cellindex(),1);            
        }

        if(init_part == 0)
        {
            if(tmp->partsize(0) == 0)
            {
                for(int j = 0; j < tmp->cellsize(); j++)
                {
                    c = celllist[tmp->cell(j)];
                    if(c->cellindex() == base->cellindex()) continue;
                    if(!c->lock())
                    {
                        int g = c->gain_num() + Pmax;
                        vector<int>::iterator it ;
                        it = find(list1[g].begin(),list1[g].end(),c->cellindex());
                        assert(it != list1[g].end());
                        list1[g].erase(it);
                        c->change_gain(-1);
                        list1[g-1].push_back(c->cellindex());
                        if(list1[g].size() == 0 && Gmax1 == g) Gmax1 = g-1; 
                    }
                }
            }
            else if(tmp->partsize(0) == 1)
            {
                c = celllist[tmp->partcellbegin(0)];
                if(!c->lock())
                { 
                    int g = c->gain_num() + Pmax;
                    vector<int>::iterator it ;
                    it = find(list0[g].begin(),list0[g].end(),c->cellindex());
                    assert(it != list0[g].end());
                    list0[g].erase(it);
                    c->change_gain(1);
                    list0[g+1].push_back(c->cellindex());
                    if (g+1 > Gmax0) Gmax0 = g+1; 
                }                
            }               
        }
        else
        {
            if(tmp->partsize(1) == 0)
            {
                for(int j = 0; j < tmp->cellsize(); j++)
                {
                    c = celllist[tmp->cell(j)];
                    if(c->cellindex() == base->cellindex()) continue;
                    if(!c->lock())
                    {
                        int g = c->gain_num() + Pmax;
                        vector<int>::iterator it ;
                        it = find(list0[g].begin(),list0[g].end(),c->cellindex());
                        assert(it != list0[g].end());
                        list0[g].erase(it);
                        c->change_gain(-1);
                        list0[g-1].push_back(c->cellindex());
                        if(list0[g].size() == 0 && Gmax0 == g) Gmax0 = g-1; 
                    }
                }
            }
            else if(tmp->partsize(1) == 1)
            {
                c = celllist[tmp->partcellbegin(0)];
                if(!c->lock())
                { 
                    int g = c->gain_num() + Pmax;
                    vector<int>::iterator it ;
                    it = find(list1[g].begin(),list1[g].end(),c->cellindex());
                    assert(it != list1[g].end());
                    list1[g].erase(it);
                    c->change_gain(1);
                    list1[g+1].push_back(c->cellindex());
                    if (g+1 > Gmax1) Gmax1 = g+1; 
                }         
            }              
        }
    }
    base->change_part();
}

void FM::iteration()
{
    int sum = 0;
    int partial_sum = 0;
    int move = 0;
    vector<int> order;
    Cell* base;
    for(int k = 0;k < nCells; k++)
    {
        if(Gmax0 >= Gmax1)
        {
            for(int i = 0; i < list0[Gmax0].size() ; i++)
            {
                base = celllist[list0[Gmax0][i]];
                if(!base->lock()) break;
            }
            sum += Gmax0-Pmax;   
        }
        else
        {
            for(int i = 0; i < list1[Gmax1].size() ; i++)
            {
                base = celllist[list1[Gmax0][i]];
                if(!base->lock()) break;
            }
            sum += Gmax1-Pmax;  
        }
        update_gain(base);
        order.push_back(base->cellindex());
        if(sum > partial_sum)
        {
            partial_sum = sum;
            move = k;
        }
    }
}

int main(int argc, char* argv[])
{
    FM* fm = new FM();
    fm->parser(argv[1]);
    fm->compute_gain();
    fm->iteration();
}
