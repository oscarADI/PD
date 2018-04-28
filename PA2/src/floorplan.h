#ifndef FLOORPLAN_H
#define FLOORPLAN_H


#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<map>
#include<ctime>
#include<cassert>
#include<cstdlib>
#include<random>
#include<cfloat>
#include"gnuplot-iostream/gnuplot-iostream.h"
#include<boost/tuple/tuple.hpp>
using namespace std;

class Block{
    public:
        Block() {}
        Block(string n,int i) {
            rotate = 0;
            name = n; 
            id = i;
            prev = next = -1;
            }
        Block(const Block& b)
        {
            rotate = b.rotate;
            name = b.name;
            id = b.id;
            width = b.width;
            height = b.height;
            next = b.next;
            prev = b.prev;
            x = b.x;
            y = b.y;
        }
        ~Block() {}

        void init(int i, double l) {(i==1) ? (width = l) : (height = l);}
        double getwidth() {return width;}
        double getheight() {return height;}
        string getname() {return name;}
        void setxy(int i, double l) {(i==0) ? (x = l) : (y = l);}
        double getxy(int i) {
            if(i==0) return x;
            else return y;
        }
        int getid() {return id;}
        void setrotate() {
            rotate ^= 1;
            swap(width,height);
        }
        int getrotate() {return rotate;}

        //for y contour
        void setnext(int b) {next = b;}
        void setprev(int b) {prev = b;}
        int nextblock() {return next;}
        int prevblock() {return prev;}
        void renew()
        {
            x = y = 0;
            next = prev = -1;
        }
        void flip() {swap(x,y);swap(width,height);}

    private:
        int id;
        double width;
        double height;
        double x;
        double y;
        string name;
        // for y contour
        int next;
        int prev;
        int rotate;

};

class Terminal{
    public:
        Terminal() {};
        Terminal(string n, int i)
        {
            name = n;
            id = i;
        }

        void init(int i, double j)
        {
            if(i == 1) x = j;
            else y = j;
        }
        string getname() {return name;}
        double getxy(int i) {return (i==0) ? x : y; }
        void flip() {swap(x,y);}
    private:
        string name;
        int id;
        double x;
        double y;
};

class Net{
    public:
        Net() {}
        Net(int i) {degree = i;}
        ~Net() {}

        void init(int i, string l) {b.push_back(i);bl.push_back(l);}
        int getdegree() {return degree;}
        int getb(int i) {return b[i];}
        string getbl(int i) {return bl[i];}

    private:
        int degree;
        vector<int> b;
        vector<string> bl;
};

class Node{
    public:
        Node() {
            data = 0;
            id = parent = left = right = -1;
        }
        Node(Block* b, int i, int p, int l, int r){
            data = b;
            id = i;
            parent = p;
            left = l;
            right = r;
        }
        Node(const Node& n)
        {
            data = n.data;
            id = n.id;
            parent = n.parent;
            left = n.left;
            right = n.right;
        }
        ~Node() {}

        Block* getdata() {return data;}
        int getinfo(int i)
        {
            if(i==0) return id;
            else if(i==1) return parent;
            else if(i==2) return left;
            else return right;
        }
        void changeinfo(int i, int j)
        {
            if(i==1) parent = j;
            else if(i==2) left = j;
            else right = j;
        }
        void renew() {
            parent = -1;
            left = -1;
            right = -1;
        }
    
    private:
        Block* data;
        int id;
        int parent;
        int left;
        int right;
};

class Result{
    public:
        Result() {}
        ~Result() {}

        void print()
        {
            cout << "Last!\n";
            for(int i = 0;i < blo.size();i++)
                cout << blo[i].getname() << " x = " << blo[i].getxy(0) << ", y = "<<blo[i].getxy(1)<<", rotate: "<<blo[i].getrotate()<<endl;
            cout << "cost = " << _cost << ", root = "<<_root<<endl<<endl;
        }
    
        vector<Block> blo;
        vector<Node> no;
        double _cost;
        double _costori;
        int _root;
        double outx;
        double outy;
        double hpwl;
};

class Floorplan{
    public:
        Floorplan() {
            alpha_base = 0;
            alpha = beta = 0;
            root = 0;
            finalX = finalY = 0;
            A_norm = W_norm = 0;
            last = 0;
            best = 0;
            cost = 0;
        }
        ~Floorplan() {}

        void seta(const char* c) {alpha = alpha_base = stod(c);}
        void setb(const char* c) {beta = stod(c);}
        void parser_block(const char*);
        void parser_net(const char*);
        void debug();
        void gnuplot();
        double HPWL();
        double Cost();
        void perturb();
        void compute_norm();
        int getRand();
        void turnback();
        void turnbackbest();

        // B*tree
        void init_tree();
        void packing(int id = 0,double xcor = 0,bool leftchild = false);
        void deletenode(int);
        void insertnode(int,int);
        void swapnode(int,int);
        void clearpack() 
        {        
            for(vector<Block*>::iterator iter = blocks.begin(); iter!=blocks.end();iter++)
                (*iter)->renew();
            finalX = 0;
            finalY = 0;
        }
        void store_result();
        void store_best();
        void SA();
        double davg(int t);
        bool fitoutline();
        void output(double, const char*);
        void check();
        
        
    
    private:
        double alpha;
        double alpha_base;
        double beta;
        double outlineX;
        double outlineY;
        int numblocks;
        int numterminals;
        int numnets;
        vector<Block*> blocks;
        vector<Net*> nets;
        vector<Terminal*> terminals;
        map<string,int> loc;
        map<string,int> ter;
        double A_norm;
        double W_norm;
        int root;
        Result* last;
        Result* best;
        double cost;

        // B*tree
        vector<Node*> bstree;
        double finalX;
        double finalY;
};

#endif