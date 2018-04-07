#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <list>
#include <time.h>
using namespace std;

class Net {
    public:
        Net() 
        {
            index = 0;
            lock0 = lock1 = 0;
        }
        Net(int n)
        {
            index = n;
            lock0 = lock1 = 0;
        }
        ~Net() {}

        // initialize
        void pushcell(int n) { cells.push_back(n); }
        int cellsize() { return cells.size(); }
        int cell(int n) { return cells[n]; }
        int netindex() { return index; }
        int partsize(int i) { if(!i) return cellsin0.size(); else return cellsin1.size(); }
        int partcellbegin(int i) { if(!i) return cellsin0[0]; else return cellsin1[0]; }
        bool lock() { return (lock0 && lock1); }
        void renew() {
            cellsin0.clear();
            cellsin1.clear();
            lock0 = 0;
            lock1 = 0;
        }
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
                // it = find(cellsin1.begin(),cellsin1.end(),n);
                // if(it != cellsin1.end()) {cout<<n<<endl;print();}
                // assert(it == cellsin1.end());
                cellsin1.push_back(n);
                lock1 = 1;
            }
            else
            {
                vector<int>::iterator it;
                //for(it = cellsin1.begin(); it!= cellsin1.end(); it++) cout<<(*it)<<" ";
                it = find(cellsin1.begin(),cellsin1.end(),n);
                assert(it != cellsin1.end());
                cellsin1.erase(it);
                // it = find(cellsin0.begin(),cellsin0.end(),n);
                // if(it != cellsin0.end()) {cout<<n<<endl;print();}
                // assert(it == cellsin0.end());
                cellsin0.push_back(n);
                lock0 = 1;
            }
            //cout<<endl<<endl;
        }

        // debug
        void print() {
            cout << "n" << index << "\n"; 
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
        int             lock0;
        int             lock1;
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
            partition = -1;
        }
        ~Cell() {}

        // initialize
        void pushnet(int n) { nets.push_back(n); }
        int Psize() { return nets.size(); }
        int cellindex() { return index; }
        void assign_part(int n) { partition = n; }
        int gain_num() { return gain; }
        void change_gain(int n) { gain += n; }
        int part_num() { return partition; }
        int net(int n) { return nets[n]; }
        int lock() {return locked;}
        void makelock() {locked = 1;}
        void unlock() {locked = gain = 0;}
        void change_part() {
            if(partition) partition = 0; 
            else partition = 1;
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
        int             locked;
        int             gain;
};

class FM {
    public:
        FM() {
            nCells = 0;
            Pmax = size0 = size1 = cutsize = 0;
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
        void            output(const char*);

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
        int                                 cutsize;
};
bool FM::balanced(int n)
{
    if(!n) 
    {
        return ((size0-1) >= lowerbound );
    }
    else
    {
        return ((size1-1) >= lowerbound );
    }
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

bool compare(Cell* i, Cell* j) {return (i->Psize()>j->Psize());}

void FM::parser (const char* argv)
{
    unsigned long start = clock();
    cout << "\nParsing" << endl;
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
                if(cell == atoi(line.c_str())) continue;
                cell = atoi(line.c_str());
                temp->pushcell(cell);
                nCells = max(nCells,cell);
            } 
        }   
    }
    infile.close();

    celllist.resize(nCells+1);
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
    lowerbound = (1-degree) * nCells / 2;
    upperbound = (1+degree) * nCells / 2;

    vector<Cell*> t = celllist;
    sort(t.begin()+1,t.end(),compare);
    
    for(vector<Cell*>::iterator iter = t.begin()+1; iter != t.end(); iter++)
    {
         Pmax = max(Pmax, (*iter)->Psize());
        if(size0 < nCells/2)
        {
            for(int i = 0;i < (*iter)->Psize();i++)
            {
                Net* n = netlist[(*iter)->net(i)];
                for(int j = 0;j < n->cellsize();j++)
                {
                    if(size0 == nCells/2) break;
                    Cell* c = celllist[n->cell(j)];
                    if(c->part_num() == -1)
                    {
                        c->assign_part(0);
                        size0++;
                    }
                }
            }
        }
        else 
        {
            if((*iter)->part_num() == -1)
            {
                (*iter)->assign_part(1);
                size1++;
            }
        }
        //cout<<"size1 = "<<size1<<endl;
        // if((*iter)->cellindex() <= nCells/2) 
        // {
        //     (*iter)->assign_part(0);
        //     size0++;
        // }
        // else 
        // {
        //     (*iter)->assign_part(1);
        //     size1++;
        // }

    }
    list0.resize(2*Pmax+1);
    list1.resize(2*Pmax+1);
    Gmax0 = Gmax1 = -Pmax;
    unsigned long end = clock();
    cout << "Number of cells : " << nCells << endl;
    cout << "Maximum Pin size : " << Pmax << endl;
    cout << "Balance : "<<lowerbound<<" <= #(G1),#(G2) <= "<<upperbound<<endl;
    cout << "Runtime of parsing : "<<(end - start)/1000000.0 <<" seconds\n"<< endl;
}

void FM::compute_gain()
{
    //cout << "compute gain !\n";
    cutsize = 0;
    for(vector<Net*>::iterator iter = netlist.begin()+1; iter != netlist.end(); iter++)
    {
        (*iter)->renew();
        int p0 = 0;
        int p1 = 0;
        for(int i = 0 ; i < (*iter)->cellsize() ; i++)
        {
            Cell* tmp = celllist[(*iter)->cell(i)];
            if(tmp->lock()) tmp->unlock();
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
        if(p0 == 0) // T(n) = 0
        {
            for(int i = 0 ; i < (*iter)->cellsize() ; i++)
            {
                Cell* tmp = celllist[(*iter)->cell(i)];
                assert(tmp->part_num()==1);
                tmp->change_gain(-1); 
            }   
        }
        else if(p1 == 0)
        {
            for(int i = 0 ; i < (*iter)->cellsize() ; i++)
            {
                Cell* tmp = celllist[(*iter)->cell(i)];
                assert(tmp->part_num()==0);
                tmp->change_gain(-1); 
            }   
        }
        if(p0 == 1)
        {
            assert((*iter)->partsize(0)==1);
            Cell* tmp = celllist[(*iter)->partcellbegin(0)];
            tmp->change_gain(1);  
        }
        if(p1 == 1)
        {
            assert((*iter)->partsize(1)==1);
            Cell* tmp = celllist[(*iter)->partcellbegin(1)];
            tmp->change_gain(1);  
        }
        if(p0 != 0 && p1 != 0) cutsize++;
    }
    size0 = size1 = 0;
    for(vector<Cell*>::iterator iter = celllist.begin()+1; iter != celllist.end(); iter++)
    {
        if((*iter)->part_num() == -1) (*iter)->assign_part(0);
        if((*iter)->part_num() == 0)
        {
            list0[(*iter)->gain_num()+Pmax].push_back((*iter)->cellindex());
            Gmax0 = max(Gmax0,(*iter)->gain_num()+Pmax);
            size0++;
        }
        else
        {
            list1[(*iter)->gain_num()+Pmax].push_back((*iter)->cellindex());
            Gmax1 = max(Gmax1,(*iter)->gain_num()+Pmax);        
            size1++;  
        }
    }
}

void FM::update_gain(Cell* base)
{
    base->makelock();
    int init_part = base->part_num();
    //vector<int> up;
    //int update[500000] = {0};
    for(int i = 0 ; i < base->Psize() ; i++)
    {
        Net* tmp = netlist[base->net(i)];
        if(tmp->lock()) 
        {
            continue;
        }
        Cell* c;
        if(init_part == 0)
        {
            if(tmp->partsize(1) == 0) // T(n) = 0
                for(int j = 0; j < tmp->cellsize(); j++)
                {
                    c = celllist[tmp->cell(j)];
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
                c = celllist[tmp->partcellbegin(1)];
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
}

void FM::iteration()
{
    unsigned long start = clock();
    cout << "Start of iteration"<<endl;
    int sum = 0;
    int partial_sum = 0;
    int move = 0;
    int de = 1;
    vector<int> order;
    Cell* base;
    while(1)
    {
        // cout << "iteration"<<de<<", ";
        move = partial_sum = sum = 0;
        //cout << "size0 = " << size0 << ", size1 = " << size1 << endl;
        for(int k = 0;k < nCells; k++)
        {
            // cout << "k = "<<k<<", Gmax0 = "<<Gmax0<<", Gmax1 = "<<Gmax1<<", Pmax = "<<Pmax<<", g = ";
            //   cout << "size0 = "<<size0<<", size1 = "<<size1<< ", lowerbound = "<<lowerbound<<endl;
            //   return;
            if( ( (Gmax0 >= Gmax1) && (balanced(0)) )  || ( (Gmax0 < Gmax1) && (!balanced(1) ) ) )
            {
                assert(list0[Gmax0].size());
                base = celllist[list0[Gmax0].front()];
                list0[Gmax0].pop_front();
                sum += Gmax0-Pmax;
                //cout<<Gmax0-Pmax<<", ";
                if(list0[Gmax0].size()==0)
                {
                    for(int i = Gmax0-1 ; ;i--)
                    {
                        if(i < 0)
                        {
                            Gmax0 = -1;
                            break;
                        }
                        if(list0[i].size())
                        {
                            Gmax0 = i;
                            break;
                        }
                    }
                }
                size0--;
                size1++;  
            }
            else
            {
                assert(list1[Gmax1].size());
                base = celllist[list1[Gmax1].front()];
                list1[Gmax1].pop_front();
                sum += Gmax1-Pmax;
                //cout<<Gmax1-Pmax<<", ";
                if(list1[Gmax1].size()==0)
                {
                    for(int i = Gmax1-1 ; ;i--)
                    {
                        if(i < 0)
                        {
                            Gmax1 = -1;
                            break;
                        }
                        if(list1[i].size())
                        {
                            Gmax1 = i;
                            break;
                        }
                    }
                }
                size0++;
                size1--;   
            }
            assert(!base->lock());
            update_gain(base);
            order.push_back(base->cellindex());
            if(sum >= partial_sum)
            {
                partial_sum = sum;
                move = k;
            }
            // cout<<"sum = "<<sum<<endl;
        }
        assert(sum == 0);
        // cout << "partial_sum = " << partial_sum << endl;
        if(partial_sum <= 0) break;
        for(int l = 0; l <= move ; l++)
        {
            Cell* c = celllist[order[l]];
            c->change_part();
        }
        compute_gain();
        order.clear();
        if(nCells > 300000 && de == 2 ) break;
        de++;
    }
    unsigned long end = clock();
    // cout << "done!"<<endl;
    cout << "End of iteration"<<endl;
    cout << "Finishing partitioning and writing output file\n";
    cout << "Runtime of iteration : "<<(end - start)/1000000.0 <<" seconds\n"<< endl;
}
void FM::output(const char* argv)
{
    unsigned long start = clock();
    ofstream outfile(argv);
    int size = 0;
    vector<int> c0,c1;
    for(vector<Cell*>::iterator iter = celllist.begin()+1 ; iter!= celllist.end() ; iter++)
    {
        if((*iter)->part_num() == 0) c0.push_back((*iter)->cellindex());
        else c1.push_back((*iter)->cellindex());
    }

    outfile << "Cutsize = " << cutsize << endl;
    outfile << "G1 "<<c0.size()<<endl;
    for(int i = 0;i < c0.size();i++)
        outfile<<"c"<<c0[i]<<" ";
    outfile<<";\n";

    outfile << "G2 "<<c1.size()<<endl;
    for(int i = 0;i < c1.size();i++)
        outfile<<"c"<<c1[i]<<" ";
    outfile<<";\n";
    
    unsigned long end = clock();
    cout << endl;
    cout << "Final Result :\n";
    cout << "   Cutsize : " << cutsize << endl;
    cout << "   Size of G1 : " << c0.size() << endl;
    cout << "   Size of G2 : " << c1.size() << endl;
    cout << "Runtime of writing output file : "<<(end - start)/1000000.0 <<" seconds\n"<< endl;
}

int main(int argc, char* argv[])
{
    unsigned long start = clock();
    FM* fm = new FM();
    fm->parser(argv[1]);
    fm->compute_gain();
    fm->iteration();
    fm->output(argv[2]);
    unsigned long end = clock();
    cout << "Total runtime : "<<(end - start)/1000000.0 <<" seconds"<< endl;
}
