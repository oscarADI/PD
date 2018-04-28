#include "floorplan.h"

void Floorplan::debug()
{
    double area = 0;
    for(int i = 0;i < numblocks;i++)
    {
        Block b = *blocks[i];
        area += b.getwidth()*b.getheight();
    }
    cout << "total area = " << area << endl;
    cout << "outline = " << outlineX*outlineY<<endl;
    cout << "ratio = " << (outlineX*outlineY-area) / (outlineX*outlineY) << endl;

    // cout << "Outline: "<<outlineX<<" "<<outlineY<<endl;
    // cout << "NumBlocks: "<<numblocks<<endl;
    // cout << "NumTerminals: "<<numterminals<<endl;
    // for(vector<Block*>::iterator iter = blocks.begin(); iter!=blocks.end(); iter++)
    //     cout<<(*iter)->getname()<<" "<<(*iter)->getwidth()<<" "<<(*iter)->getheight()<<endl;
    // for(vector<Terminal*>::iterator iter = terminals.begin();iter!=terminals.end();iter++)
    //     cout<<(*iter)->getname()<<" terminal "<<(*iter)->getxy(0)<<" "<<(*iter)->getxy(1)<<endl;
    // cout << "NumNets: " << numnets << endl;
    //     for(vector<Net*>::iterator iter = nets.begin(); iter!=nets.end(); iter++)
    //     {
    //         cout<<"NetDegree: "<<(*iter)->getdegree()<<"\n";
    //         for(int i = 0;i < (*iter)->getdegree(); i++)
    //             cout<<(*iter)->getbl(i)<<endl;
    //     }
    // store_result();
    // for(int i = 0;i < numblocks;i++)
    // {
    //     Block b = best->blo[i];
    //     cout << b.getname() << " x = " << b.getxy(0) << ", y = "<<b.getxy(1)<<", rotate: "<<b.getrotate()<<endl;
    // }
    // cout << "cost = " << cost << endl;
    // cout << endl;
    // last->print();

    // for(int j = 0;j < 5;j++)
    // {
    //     cout << "j = " << j << endl;
    //     perturb();
    //     cout << "Before!\n";
    //     for(int i = 0;i < numblocks;i++)
    //     {
    //         Block b = *blocks[i];
    //         cout << b.getname() << " x = " << b.getxy(0) << ", y = "<<b.getxy(1)<<", rotate: "<<b.getrotate()<<endl;
    //     }
    //     cout << "cost = " << cost << endl;
    //     cout << endl;
    //     if(j%2 == 0)
    //     {
    //         store_result();
    //         last->print();
    //     }
    //     else
    //     {
    //         turnback();
    //     }
    //     cout << "After!\n";
    //     for(int i = 0;i < numblocks;i++)
    //     {
    //         Block b = *blocks[i];
    //         cout << b.getname() << " x = " << b.getxy(0) << ", y = "<<b.getxy(1)<<", rotate: "<<b.getrotate()<<endl;
    //     }
    //     cout << "cost = " << cost << endl;
    //     cout << endl;
    // }
}

void Floorplan::parser_block(const char* argv)
{
    cout << "Parsing block file\n";
    ifstream infile;
    infile.open(argv,ifstream::in);
    if(!infile.is_open())
    {
        cerr << "Input file doesn't exist !" << endl;
        return;
    }

    string line;
    while(infile >> line){
        if(line == "Outline:")
        {
            infile >> line;
            outlineX = stod(line);
            infile >> line;
            outlineY = stod(line);
        }
        else if(line == "NumBlocks:")
        {
            infile >> line;
            numblocks = atoi(line.c_str());
        }
        else if(line == "NumTerminals:")
        {
            infile >> line;
            numterminals = atoi(line.c_str());
            break;
        }
    }
    for(int i = 0;i < numblocks;i++)
    {
        infile >> line;
        Block* temp = new Block(line,blocks.size());
        infile >> line;
        temp->init(1,stod(line));
        infile >> line;
        temp->init(2,stod(line));
        loc[temp->getname()] = blocks.size();
        blocks.push_back(temp);
    }
    for(int i = 0;i < numterminals;i++)
    {
        infile >> line;
        Terminal* temp = new Terminal(line, terminals.size());
        infile >> line;
        infile >> line;
        temp->init(1,stod(line));
        infile >> line;
        temp->init(2,stod(line));
        ter[temp->getname()] = terminals.size();
        terminals.push_back(temp);
    }
}

void Floorplan::parser_net(const char* argv)
{
    cout << "Parsing net file\n";
    ifstream infile;
    infile.open(argv,ifstream::in);
    if(!infile.is_open())
    {
        cerr << "Input file doesn't exist !" << endl;
        return;
    }

    string line;
    while(infile >> line)
    {
        if(line == "NumNets:")
        {
            infile >> line;
            numnets = atoi(line.c_str());
        }
        else if(line == "NetDegree:")
        {
            infile >> line;
            int d = atoi(line.c_str());
            Net* temp = new Net(d);
            for(int i = 0;i < d;i++)
            {
                infile >> line;
                if(loc.find(line) != loc.end())
                    temp->init(loc[line],line);
                else if(ter.find(line) != loc.end())
                    temp->init(ter[line],line);
                else 
                {
                    cout << line << endl;
                    assert(false);
                }
            }
            nets.push_back(temp);
        }
    }
}

void Floorplan::init_tree()
{
    for(int i = 0;i < numblocks; i++)
    {
        int p,l,r;
        p = (i-1)/2;
        l = 2*i+1;
        r = 2*i+2;
        if(i==0) p = -1;
        if(l>=numblocks) l = -1;
        if(r>=numblocks) r = -1;
        Node* n = new Node(blocks[i],i,p,l,r);
        bstree.push_back(n);
    }
}

void Floorplan::packing(int id, double xcor, bool leftchild)
{
    // cout << "packing";
    // cout << " ID = "<<id;
    assert(id != -1);
    Node* n = bstree[id];
    Block* b = n->getdata();
    b->setxy(0,xcor); // Set b->x = xcor
    if(id == root)
    {
        b->setxy(0,0);
        b->setxy(1,0);
        finalX = 0;
        finalY = 0;
        //cout << ", x = 0, y = 0, width = "<<b->getwidth() << ", height = "<<b->getheight()<<endl;
        if(n->getinfo(2) != -1) packing(n->getinfo(2),b->getwidth(),1); //pack left child
        if(n->getinfo(3) != -1) packing(n->getinfo(3),0,0); // pack right child
        return;
    }
    assert(n->getinfo(1) != -1);
    Block* parent = bstree[n->getinfo(1)]->getdata();
    // cout << " parent = "<<parent->getid() << endl;

    //for y
    int ref;
    if(leftchild) ref = parent->nextblock(); // reference : parent->next
    else ref = parent->getid(); // reference : parent
    double ymax = 0;
    for(int s = ref; s != -1; s = blocks[s]->nextblock())
    {
        Block* start = blocks[s];
        double leftx = start->getxy(0);
        double rightx = leftx + start->getwidth();

        if(xcor < rightx && xcor+b->getwidth() > rightx)
        {
            ymax = max(ymax, start->getxy(1)+start->getheight());
        }
        else if(xcor < rightx && xcor+b->getwidth() <= rightx)
        {
            ymax = max(ymax, start->getxy(1)+start->getheight());
            if(xcor+b->getwidth() == rightx) b->setnext(blocks[s]->nextblock());
            else b->setnext(s);
            break;
        }
        // else
        // {
        //     b->setnext(s);
        //     break;
        // }
    }
    assert(b->getid()!=-1);
    assert(parent->getid()!=-1);
    // if(xcor == 0)
    // {   
    //     //parent->setnext(b->getid());
    //     b->setprev(parent->getid());
    //     //b->setnext(-1);
    // }
    if(leftchild) // left child
    {
        b->setprev(parent->getid());
        // cout << "parent = "<<parent->getid() << endl;
        // cout << "xcor1 = " << xcor << endl;
        assert(parent->getid() != -1);
        parent->setnext(b->getid());
    }
    else// right child
    {
        b->setprev(parent->prevblock());
        // cout << parent->prevblock() << endl;
        // cout << "xcor2 = " << xcor << endl;
        if(xcor != 0)
        {
            assert(parent->prevblock()!=-1);
            blocks[parent->prevblock()]->setnext(b->getid());
        }
    }  
    b->setxy(1,ymax);
    finalX = max(finalX, xcor+b->getwidth());
    finalY = max(finalY,ymax+b->getheight());
    //cout << ", x = " << b->getxy(0) << ", y = " << b->getxy(1) <<", width = "<<b->getwidth()<<", height = "<<b->getheight()<< endl;
    if(n->getinfo(2) != -1) packing(n->getinfo(2),xcor+b->getwidth(),1); //pack left child
    if(n->getinfo(3) != -1) packing(n->getinfo(3),xcor,0); // pack right child
}

double Floorplan::HPWL()
{
    double hpwl = 0;
    for(int i = 0;i < numnets; i++)
    {
        Net* n = nets[i];
        double minX = outlineX, minY = outlineY;
        double maxX = 0, maxY = 0;
        for(int j = 0;j < n->getdegree();j++)
        {
            if(loc.find(n->getbl(j))!=loc.end())
            {
                Block* b = blocks[n->getb(j)];
                minX = min(minX,b->getxy(0)+b->getwidth()/2);
                minY = min(minY,b->getxy(1)+b->getheight()/2);
                maxX = max(maxX,b->getxy(0)+b->getwidth()/2);
                maxY = max(maxY,b->getxy(1)+b->getheight()/2);
            }
            else
            {
                Terminal* t = terminals[n->getb(j)];
                minX = min(minX,t->getxy(0));
                minY = min(minY,t->getxy(1));
                maxX = max(maxX,t->getxy(0));
                maxY = max(maxY,t->getxy(1));                
            }
        }
        // cout << minX << endl;
        // cout << minY << endl;
        // cout << maxX << endl;
        // cout << maxY << endl;
        // cout << "hpwl for net"<<i<<" = "<< (maxX-minX)+(maxY-minY)<<endl;
        hpwl += (maxX-minX)+(maxY-minY);
    }
    // cout << "W = " << hpwl << endl;
    // cout << "X = "<<finalX<<", Y = "<<finalY<<", A = " << finalX*finalY << endl;
    return hpwl;
}

double Floorplan::Cost()
{
    double Wr = outlineX;
    double Hr = outlineY;
    double Wrs = finalX;
    double Hrs = finalY;
    if(Wr < Hr) swap(Wr,Hr);
    if(Wrs < Hrs) swap(Wrs,Hrs);
    double R = Hrs/Wrs - Hr/Wr;

    cost =  (alpha*finalX*finalY/A_norm) + (beta*HPWL()/W_norm) + 
            (1-alpha-beta)*R*R;
    return cost;
}
int Floorplan::getRand()
{
    random_device rd;
    uniform_int_distribution<int> dist(0, numblocks - 1);
    return dist(rd);
}
void Floorplan::perturb()
{
    // int op = 1; 
    int op = getRand()%3;
    if(op == 0)
    {
        
        int target = getRand()%numblocks;
        //cout << "rotate "<<target<<" ";
        blocks[target]->setrotate();
    }
    else if(op == 1)
    {
        int t1 = getRand()%numblocks;
        int t2 = getRand()%numblocks;
        while(t1 == t2) t2 = getRand()%numblocks;
        //cout << "delete "<<t1<< " to " << t2 <<" ";
        deletenode(t1);
        insertnode(t1,t2);   
    }
    else
    {
        int t1 = getRand()%numblocks;
        int t2 = getRand()%numblocks;
        while(t1 == t2) t2 = getRand()%numblocks;
        swapnode(t1,t2);
        //cout << "swap "<<t1<<"&"<<t2<<" ";
    }
    clearpack();
    packing(root);
}

void Floorplan::deletenode(int id)
{
    Node* n = bstree[id];
    int p = n->getinfo(1); // parent
    int l = n->getinfo(2); // left child
    int r = n->getinfo(3); // right child
    if(l == -1 && r == -1)
    {
        if(bstree[p]->getinfo(2) == id) bstree[p]->changeinfo(2,-1);
        else if(bstree[p]->getinfo(3) == id) bstree[p]->changeinfo(3,-1);
        else 
        {
            cout << bstree[p]->getinfo(2) << endl;
            cout << bstree[p]->getinfo(3) << endl;
            cerr<<"wrong delete\n";
        }
        n->renew();
    }
    else if(l == -1 && r != -1)
    {
        if(p != -1)
        {
            if(bstree[p]->getinfo(2) == id) bstree[p]->changeinfo(2,r);
            else if(bstree[p]->getinfo(3) == id) bstree[p]->changeinfo(3,r);
            else 
            {
                cout << bstree[p]->getinfo(2) << endl;
                cout << bstree[p]->getinfo(3) << endl;
                cerr<<"wrong delete\n";
            }
        }
        else root = r;
        bstree[r]->changeinfo(1,p);
        n->renew();
    }
    else if(l != -1 && r == -1)
    {
        if(p != -1)
        {
            if(bstree[p]->getinfo(2) == id) bstree[p]->changeinfo(2,l);
            else if(bstree[p]->getinfo(3) == id) bstree[p]->changeinfo(3,l);
            else 
            {
                cout << bstree[p]->getinfo(2) << endl;
                cout << bstree[p]->getinfo(3) << endl;
                cerr<<"wrong delete\n";
            }
        }
        else root = l;
        bstree[l]->changeinfo(1,p);
        n->renew();
    }
    else
    {
        int s = rand() % 2;
        if(s == 0) swapnode(id,l);
        else swapnode(id,r);
        deletenode(id);
    }
}

void Floorplan::insertnode(int id, int parent)
{
    // cout << "insert "<<id<<" to "<<parent<<endl;
    Node* child = bstree[id];
    Node* p = bstree[parent];
    child->changeinfo(1,parent);

    int s = rand() % 2;
    int c = -1;
    if(s == 0) //left of parent
    {
        c = p->getinfo(2);
        if(c != -1)
        {
            bstree[c]->changeinfo(1,id);
        }
        p->changeinfo(2,id);
        child->changeinfo(2,c);
    }
    else
    {
        c = p->getinfo(3);
        if(c != -1)
        {
            bstree[c]->changeinfo(1,id);
        }
        p->changeinfo(3,id);
        child->changeinfo(3,c);
    }

    // int rnd = rand()%2;
    // if(rnd == 0) child->changeinfo(2,c);
    // else child->changeinfo(3,c);
}

void Floorplan::swapnode(int n1, int n2)
{
    // cout << "swapnode "<<n1<<"&"<<n2<<endl;
    Node* nd1 = bstree[n1];
    Node* nd2 = bstree[n2];
    int p1 = nd1->getinfo(1);
    int l1 = nd1->getinfo(2);
    int r1 = nd1->getinfo(3);
    int p2 = nd2->getinfo(1);
    int l2 = nd2->getinfo(2);
    int r2 = nd2->getinfo(3);
    // n2->n1
    nd2->changeinfo(1,p1);
    if(n1 != root)
    {
        assert(p1 != -1);
        if(bstree[p1]->getinfo(2) == n1) bstree[p1]->changeinfo(2,n2);
        else bstree[p1]->changeinfo(3,n2);
    }
    nd2->changeinfo(2,l1);
    if(l1 != -1) bstree[l1]->changeinfo(1,n2);
    nd2->changeinfo(3,r1);
    if(r1 != -1) bstree[r1]->changeinfo(1,n2);

    // n1->n2
    nd1->changeinfo(1,p2);
    if(n2 != root)
    {
        assert(p2 != -1);
        if(bstree[p2]->getinfo(2) == n2) bstree[p2]->changeinfo(2,n1);
        else bstree[p2]->changeinfo(3,n1);
    }
    nd1->changeinfo(2,l2);
    if(l2 != -1) bstree[l2]->changeinfo(1,n1);
    nd1->changeinfo(3,r2);
    if(r2 != -1) bstree[r2]->changeinfo(1,n1);    
    if(p1 == n2) // n2 = parent(n1)
    {
        nd2->changeinfo(1,n1);
        nd1->changeinfo(1,p2);
        if(l2 == n1) nd1->changeinfo(2,n2);
        else if(r2 == n1) nd1->changeinfo(3,n2);
        else
        {
            cerr << "wrong swap\n";
        }
    }
    else if(p2 == n1)
    {
        nd1->changeinfo(1,n2);
        nd2->changeinfo(1,p1);
        if(l1 == n2) nd2->changeinfo(2,n1);
        else if(r1 == n2) nd2->changeinfo(3,n1);
        else
        {
            cerr << "wrong swap2\n";
        }
    }
    if(n1 == root) root = n2;
    else if(n2 == root) root = n1;
}

void Floorplan::compute_norm()
{
    int t = 20 * numblocks;
    double w = 0;
    double a = 0;
    for(int i = 0;i < t;i++)
    {
        perturb();
        w += HPWL();
        a += finalX*finalY;
        if(fitoutline()) cout << "Fit!!\n";
    }
    W_norm = w/t;
    A_norm = a/t;
    cout << "W_norm = "<<W_norm<<", A_norm = "<<A_norm<<endl;
}

void Floorplan::store_result()
{
    // cout << "storeResult\n";
    last->blo.clear();
    last->blo.resize(numblocks);
    last->no.clear();
    last->no.resize(numblocks);
    for(int i = 0;i < numblocks;i++)
    {
        last->blo[i] = *blocks[i];
        last->no[i] = *bstree[i];
    }
    last->_cost = cost;
    last->_root = root;
    last->outx = finalX;
    last->outy = finalY;
}

void Floorplan::store_best()
{
    best->blo.clear();
    best->blo.resize(numblocks);
    best->no.clear();
    best->no.resize(numblocks);
    for(int i = 0;i < numblocks;i++)
    {
        best->blo[i] = *blocks[i];
        best->no[i] = *bstree[i];
    }
    best->_cost = cost;
    best->_root = root;
    best->_costori = (alpha_base*finalX*finalY/A_norm) + (1-alpha_base)*HPWL()/W_norm;
    best->outx = finalX;
    best->outy = finalY;
    best->hpwl = HPWL();
}
void Floorplan::turnback()
{
    // cout << "turnback\n";
    for(int i = 0;i < numblocks;i++)
    {
        *blocks[i] = last->blo[i];
        *bstree[i] = last->no[i];
    }
    root = last->_root;
    cost = last->_cost;
    // packing(root);
}
void Floorplan::turnbackbest()
{
    for(int i = 0;i < numblocks;i++)
    {
        *blocks[i] = best->blo[i];
        *bstree[i] = best->no[i];
    }
    root = best->_root;
    clearpack();
    packing(root);
}
void Floorplan::SA()
{
    cout << "SA\n";
    // init
    double p = 0.99;
    int k = 7;
    double ep = 0.00001;
    double reject = 0;
    int N = k*numblocks;
    double MT = 0,uphill = 0;
    Cost();
    last = new Result();
    best = new Result();
    store_result();
    store_best();
    best->_cost = DBL_MAX;
    double T1 = fabs(davg(N)/log(p));
    double T = T1;
    int n = 1;

    // cost components
    alpha = alpha_base;
    beta = 0;

    int change = 0;
    int feasible = 0;
    double count = 0;
    double avgcost = 0;
  
    //simulate
    cout << "base = "<<alpha_base<<endl;
    //return;
    while(1)
    {
        MT = uphill = reject = 0;
        n++;
        avgcost = count = 0;
        // cout << "n = " << n ;
        while(1)
        {
            double prev = cost;
            // cout << "prev = "<<prev;
            perturb();
            Cost();
            // cout << "  best : "<<best->_cost << "\t cost = "<<cost<<endl;
            MT++;
            count++;
            double dcost = cost - prev;
            avgcost += abs(dcost);
            p = min(1.0,exp(-1*dcost/T));
            if(dcost <= 0 || (getRand()%100+1)/100 < p)
            {
                if(dcost > 0) uphill++;
                store_result();
                if(cost < best->_cost) 
                {
                    if(fitoutline())
                    {
                        cout << "new best\n";
                        store_best();
                        change++;
                        alpha = alpha_base + (1-alpha_base)*change/2/N;
                    }
                }
            }
            else 
            {
                turnback();
                reject++;
            }
            if(uphill > N || MT > 2*N) break;
        }
        // update T
        double ac = avgcost / count;
        // cout << " ac = "<<ac <<endl;
        if(n >= 2 && n <= k) T = T1*ac/n/100;
        else T = T1*ac/n;

        if(reject/MT > 0.95 || T < ep ){ 
            if(change != 0) break;
            n = 1;
            T = T1;
            p = 0.99;
            alpha = alpha_base;
            feasible++;
            cout << "feasible = " << feasible << endl;   
        }
    }
    cout << "Best!\n";
    for(int i = 0;i < numblocks;i++)
    {
        Block b = best->blo[i];
        cout << b.getname() << " x = " << b.getxy(0) << ", y = "<<b.getxy(1)<<", rotate: "<<b.getrotate()<<endl;
    }
    cout << "cost = " << best->_costori << endl;
    cout << "area = " << best->outx*best->outy << endl;
    cout << "wirelength = " << best->hpwl << endl;
    cout << endl;
    turnbackbest();
    for(int i = 0;i < numblocks;i++)
    {
        Block b = *blocks[i];
        cout << b.getname() << " x = " << b.getxy(0) << ", y = "<<b.getxy(1)<<", rotate: "<<b.getrotate()<<endl;
    }
    cout << "cost = " << Cost() << endl;
    cout << "area = " << finalX*finalY << endl;
    cout << "wirelength = " << HPWL() << endl;
    cout << endl;

}
double Floorplan::davg(int t)
{
    double c = 0, uphill = 0;
    double u = 0;
    double prev = cost;

    while(uphill == 0)
    {
        u = 0;
        for(int i = 0;i < t;i++)
        {
            perturb();
            //packing();
            Cost();
            if(cost - prev > 0)
            {
                uphill += cost-prev;
                prev = cost;
                u++;
            }
        }
    }
    turnback();
    return uphill/u;
}
bool Floorplan::fitoutline()
{
    // cout << outlineX << " " << outlineY << " " << finalX << " " << finalY << endl;
    double x1 = max(outlineX,outlineY);
    double y1 = min(outlineX,outlineY);
    double x2 = max(finalX,finalY);
    double y2 = min(finalX,finalY);

    return (x2 <= x1) && (y2 <= y1);
}
// void Floorplan::gnuplot()
// {
//     Gnuplot plot;
//     double xr = 1.2*max(outlineX,outlineY), yr = xr;
//     plot << "set xrange [" << -0.1*xr << ":" << xr << "]" << endl;
//     plot << "set yrange [" << -0.1*xr << ":" << xr << "]" << endl;
//     plot << "set size ratio -1" << endl;
//     plot << "set object 2 rect from 0,0 to " << outlineX
//          << "," << outlineY << "fc rgb \"#CCFFFF\" back" << endl;
//     for (int i = 0; i < numblocks; ++i) {
//         Block b = best->blo[i];
//         plot << " set label \"" << b.getname() << "\" at "
//              << b.getxy(0)+b.getwidth()/2 << "," << b.getxy(1)+b.getheight()/2 << endl;
//         plot << " set object " << i + 3 << " rect from "
//              << b.getxy(0) << ","
//              << b.getxy(1) << " to "
//              << b.getxy(0)+b.getwidth() << ","
//              << b.getxy(1)+b.getheight() << " fc rgb \"green\" " << endl;
//     }
//     for (int i = 0; i < numterminals; ++i) {
//         Terminal* t = terminals[i];
//         plot << " set object " << i + 3 + numblocks << " circle at "
//              << t->getxy(0) << "," << t->getxy(1)
//              << " size scr 0.01 fc rgb \"navy\"" << endl;
//     }
//     plot << " plot \'-\' using 1:2 t \"\" with line" << endl;
//     plot << " 0 0" << endl;
//     plot << " e" << endl;
//     plot << " pause -1" << endl;
// }
void Floorplan::output(double t, const char* argv)
{
    ofstream outfile(argv);
    // bool f = false;
    // if(best->outx > outlineX || best->outy > outlineY) 
    // {
    //     cout << best->outx << " " << best->outy<<endl;
    //     cout << outlineX << " " << outlineY << endl;
    //     f = true;
    // }
    if(best->outx > outlineX || best->outy > outlineY)
    {
        // swap(best->outx,best->outy);
        for (int i = 0; i < numblocks; ++i) 
        {
            best->blo[i].flip();
            if(i < numterminals) terminals[i]->flip();
        }
    }

    // <final cost>
    outfile << best->_cost << endl;
    // <total wirelength>
    outfile << best->hpwl << endl;
    // <chip_area>
    outfile << best->outx*best->outy << endl;
    // <chip_width><chip_height>
    outfile << best->outx << " " << best->outy << endl;
    // else outfile << best->outy << " " << best->outx << endl;
    // <program_runtime>
    outfile << t << endl;
    double x = 0,y = 0;
    for(int i = 0;i < numblocks;i++)
    {
        Block b = best->blo[i];
        outfile << b.getname()<<"\t";
        if(1)
        {
            outfile << b.getxy(0) << "\t" << b.getxy(1) << "\t";
            outfile << b.getxy(0)+b.getwidth() << "\t" << b.getxy(1)+b.getheight() << endl;
            y = max(y, b.getxy(1)+b.getheight());
            x = max(x, b.getxy(0)+b.getwidth()); 
        }
        // else
        // {
        //     outfile << b.getxy(1) << "\t" << b.getxy(0) << "\t";
        //     outfile << b.getxy(1)+b.getheight() << "\t" << b.getxy(0)+b.getwidth() << endl; 
        //     x = max(x, b.getxy(1)+b.getheight());
        //     y = max(y, b.getxy(0)+b.getwidth()); 
        // }    
    }
    if (x != best->outx || y != best->outy)
    {
        cout << "x = "<<x << " best->outx = " << best->outx<<endl;
        cout << "y = "<<y << " best->outy = " << best->outy<<endl;
    }
}

void Floorplan::check()
{
    for(int i = 0;i < numblocks-1;i++)
    {
        for (int j = 0;j < numblocks;j++)
        {
            double diffx = best->blo[i].getxy(0) - best->blo[j].getxy(0);
            double diffy = best->blo[i].getxy(1) - best->blo[j].getxy(1);
            if(diffx >= 0 && diffy >= 0)
            {
                if(diffx < best->blo[j].getwidth() || diffy < best->blo[j].getheight())
                {
                    cerr << "wrong" << endl;
                }
            }
            else if (diffx >= 0 && diffy < 0)
            {
                if(diffx < best->blo[j].getwidth() || abs(diffy) < best->blo[i].getheight())
                {
                    cerr << "wrong" << endl;
                }
            }
            else if (diffx < 0 && diffy >= 0)
            {
                if(abs(diffx) < best->blo[i].getwidth() || diffy < best->blo[j].getheight())
                {
                    cerr << "wrong" << endl;
                }
            }
            else
            {
                if(abs(diffx) < best->blo[i].getwidth() || abs(diffy) < best->blo[i].getheight())
                {
                    cerr << "wrong" << endl;
                }                
            }
        }
    }
}