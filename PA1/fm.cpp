#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <list>
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
            //cout<<"n = "<<n<<endl;
            if(!i)
            {
                vector<int>::iterator it;
                //for(it = cellsin0.begin(); it!= cellsin0.end(); it++) cout<<(*it)<<" ";
                it = find(cellsin0.begin(),cellsin0.end(),n);
                if(it == cellsin0.end()){
                    cout << "c" <<n<<endl;
                    print();
                }
                assert(it != cellsin0.end());
                cellsin0.erase(it);
                cellsin1.push_back(n);
            }
            else
            {
                vector<int>::iterator it;
                //for(it = cellsin1.begin(); it!= cellsin1.end(); it++) cout<<(*it)<<" ";
                it = find(cellsin1.begin(),cellsin1.end(),n);
                assert(it != cellsin1.end());
                cellsin1.erase(it);
                cellsin0.push_back(n);
            }
            //cout<<endl<<endl;
        }

        // debug
        void print() {
            cout << "n" << index << " "; 
            for(vector<int>::iterator iter = cellsin0.begin(); iter != cellsin0.end(); iter++)
                cout << "c" << *iter << " ";
            cout<<endl;
            for(vector<int>::iterator iter = cellsin1.begin(); iter != cellsin1.end(); iter++)
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
            Pmax = size0 = size1 = 0;
            degree = lowerbound = upperbound = 0.0;
            Net* net = new Net(0);
            netlist.push_back(net);
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
        void            print_nets() {
            for(vector<Net*>::iterator iter = netlist.begin(); iter != netlist.end(); iter++)
            {
                (*iter)->print();
                if((*iter)->netindex() == 54293) break;
            }
        }

    private:
        int                                 nCells;
        vector<Net*>                        netlist;
        vector<Cell*>                       celllist;
        double                              degree;
        vector< list<int> >                 list0;
        vector< list<int> >                 list1;
        int                                 Pmax;
        double                              lowerbound;
        double                              upperbound;
        int                                 Gmax0;
        int                                 Gmax1;
        int                                 size0;
        int                                 size1;
};
bool FM::balanced(int n)
{
    int s0 = size0;
    int s1 = size1;
    if(n) 
    {
        s0+=1;
        s1-=1;
    }
    else
    {
        s0-=1;
        s1+=1;
    }
    return ((s0>=lowerbound)&&(s0<=upperbound)&&(s1>=lowerbound)&&(s1<=upperbound));
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
    cout << "parser !" << endl;
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
            current_net++;
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
            //if((*iter)->netindex() == 0) cout<<"0!\n";
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
                size0++;
            }
            else 
            {
                (*iter)->assign_part(1);
                size1++;
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
    cout << "conpute gain !\n";
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
            Cell* tmp = celllist[(*iter)->partcellbegin(0)];
            tmp->change_gain(1);  
        }
        if(p1 == 1)
        {
            Cell* tmp = celllist[(*iter)->partcellbegin(1)];
            tmp->change_gain(1);  
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
    //cout << "update gain !\n";
    int init_part = base->part_num();
    //cout << "init_part = " << init_part << endl;
    for(int i = 0 ; i < base->Psize() ; i++)
    {
        Net* tmp = netlist[base->net(i)];
        //cout << "net = "<<tmp->netindex()<<endl;
        Cell* c;
        if(init_part == 0)
        {
            //cout << "2\n";
            if(tmp->partsize(1) == 0) // T(n) = 0
                for(int j = 0; j < tmp->cellsize(); j++)
                {
                    c = celllist[tmp->cell(j)];
                    if(c->cellindex() == base->cellindex()) continue;
                    if(!c->lock())
                    {
                        int g = c->gain_num() + Pmax;
                        list0[g].remove(c->cellindex());
                        c->change_gain(1);
                        list0[g+1].push_back(c->cellindex());
                        if (g+1 > Gmax0) Gmax0 = g+1;
                    }

                }
            else if(tmp->partsize(1) == 1) // T(n) = 1
            {
                c = celllist[tmp->partcellbegin(1)];
                if(!c->lock())
                { 
                    int g = c->gain_num() + Pmax;
                    list1[g].remove(c->cellindex());
                    c->change_gain(-1);
                    list1[g-1].push_back(c->cellindex());
                    if(list1[g].size() == 0 && Gmax1 == g) Gmax1 = g-1; 
                }
            }
            tmp->movepart(base->cellindex(),0);
        }
        else  //init part = 1
        {
            if(tmp->partsize(0) == 0) // T(n) = 0
            {
                for(int j = 0; j < tmp->cellsize(); j++)
                {
                    c = celllist[tmp->cell(j)];
                    if(c->cellindex() == base->cellindex()) continue;
                    if(!c->lock())
                    {
                        int g = c->gain_num() + Pmax;
                        list1[g].remove(c->cellindex());
                        c->change_gain(1);
                        list1[g+1].push_back(c->cellindex());
                        if (g+1 > Gmax1) Gmax1 = g+1;
                    }
                }
            }
            else if(tmp->partsize(0) == 1) // T(n) = 1
            {
                c = celllist[tmp->partcellbegin(0)];
                if(!c->lock())
                { 
                    int g = c->gain_num() + Pmax;
                    list0[g].remove(c->cellindex());
                    c->change_gain(-1);
                    list0[g-1].push_back(c->cellindex());
                    if(list0[g].size() == 0 && Gmax0 == g) Gmax0 = g-1; 
                }                
            }
            tmp->movepart(base->cellindex(),1);            
        }

        if(init_part == 0)
        {
            if(tmp->partsize(0) == 0) // F(n) = 0
            {
                //cout << "1\n";
                for(int j = 0; j < tmp->cellsize(); j++)
                {
                    c = celllist[tmp->cell(j)];
                    if(c->cellindex() == base->cellindex()) continue;
                    if(!c->lock())
                    {
                        int g = c->gain_num() + Pmax;
                        list1[g].remove(c->cellindex());
                        c->change_gain(-1);
                        list1[g-1].push_back(c->cellindex());
                        if(list1[g].size() == 0 && Gmax1 == g) Gmax1 = g-1; 
                    }
                }
            }
            else if(tmp->partsize(0) == 1)
            {
                //cout << "11\n";
                c = celllist[tmp->partcellbegin(0)];
                if(!c->lock())
                { 
                    int g = c->gain_num() + Pmax;
                    list0[g].remove(c->cellindex());
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
                //cout << "111\n";
                for(int j = 0; j < tmp->cellsize(); j++)
                {
                    c = celllist[tmp->cell(j)];
                    if(c->cellindex() == base->cellindex()) continue;
                    if(!c->lock())
                    {
                        int g = c->gain_num() + Pmax;
                        list0[g].remove(c->cellindex());
                        c->change_gain(-1);
                        list0[g-1].push_back(c->cellindex());
                        if(list0[g].size() == 0 && Gmax0 == g) Gmax0 = g-1; 
                    }
                }
            }
            else if(tmp->partsize(1) == 1)
            {
                //cout << "1111\n";
                c = celllist[tmp->partcellbegin(0)];
                if(!c->lock())
                { 
                    int g = c->gain_num() + Pmax;
                    list1[g].remove(c->cellindex());
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
    cout << "iteration!\n";
    int sum = 0;
    int partial_sum = 0;
    int move = 0;
    vector<int> order;
    Cell* base;
    for(int k = 0;k < nCells; k++)
    {
        cout << "k = "<<k<<", ";
        if((Gmax0 >= Gmax1 && balanced(0)) || (Gmax0 < Gmax1 && !balanced(1)))
        {
            base = celllist[list0[Gmax0].front()];
            list0[Gmax0].pop_front();
            if(list0[Gmax0].size()==0)
            {
                for(int i = Gmax0-1 ; i > -Pmax ;i--)
                {
                    if(list0[i].size())
                    {
                        Gmax0 = i;
                        break;
                    }
                }
            }
            sum += Gmax0-Pmax;
            size0--;
            size1++;   
        }
        else
        {
            base = celllist[list1[Gmax1].front()];
            list1[Gmax1].pop_front();
            if(list1[Gmax1].size()==0)
            {
                for(int i = Gmax1-1 ; i > -Pmax ;i--)
                {
                    if(list1[i].size())
                    {
                        Gmax1 = i;
                        break;
                    }
                }
            }
            sum += Gmax1-Pmax;
            size0++;
            size1--;  
        }
        update_gain(base);
        order.push_back(base->cellindex());
        if(sum > partial_sum)
        {
            partial_sum = sum;
            move = k;
        }
        cout<<"sum = "<<sum<<endl;
    }
}

int main(int argc, char* argv[])
{
    FM* fm = new FM();
    fm->parser(argv[1]);
    fm->compute_gain();
    fm->iteration();
}
