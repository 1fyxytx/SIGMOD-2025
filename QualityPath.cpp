#include <omp.h>
#include <time.h>
#include <iostream>
#include <stack>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <utility>
#include <set>
#include <vector>
#include <queue>
#include <map>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <bitset>
#include <math.h>

#define MAXINT ((unsigned) 4294967295)

using namespace std;

struct trielem{
    unsigned vid, wval, flgv;

    trielem(){
        vid=0, wval=0, flgv=0;
    }

    trielem(unsigned val1, unsigned val2, unsigned val3){
        vid=val1, wval=val2, flgv=val3;
    }
};
    
bool cmp(trielem e1, trielem e2){
    if (e1.wval > e2.wval)
        return true;
    else if (e1.wval == e2.wval){
        if (e1.flgv > e2.flgv)
            return true;
        else 
            return false;
    }else
        return false;
}


struct subt{
    long long cnts;
    vector<unsigned> que, lab;
    
    subt(){
        cnts = 0;
    }

    subt(long long val, vector<unsigned>& val2, 
            vector<unsigned>& val3){
        cnts = val;
        que.insert(que.end(), val2.begin(), val2.end());
        lab.insert(lab.end(), val3.begin(), val3.end());
    }
};

bool cmp1(subt e1, subt e2){
    if (e1.cnts >= e2.cnts)
        return true;
    else
        return false;
}

vector<vector<unsigned> > con, label, labelR, Labs; // id+weight
vector<vector<trielem> > conR; //  d+d_s
vector<vector< vector<long long> > > SearchPath;

vector<vector<subt> > TaskList;
vector<subt> subtasks;

unsigned MAXDIS, MAXMOV, MASK;
int totalV = 0, weig = 20, threads = 20;
int  src, dst, aplace, bplace, fflg = 0;
unsigned CurW = 0, dmax = 0, tag=0; // CurW 表示当前搜索情况下的最大下限值

unordered_map<unsigned, unsigned> d2w;
vector<unsigned> w2d;

long long paths = 0, uness = 0, UpperBound = 0, Wavg = 0, MaxLoad;

vector<long long> PathCnt;
vector<pair<int, int> > v2degree;


vector<pair<int, int> > QPairs;
vector<pair<int, float> > D2T;



void InitialP(){
    vector<vector<unsigned> >().swap(label),
    vector<vector<unsigned> >().swap(labelR), 
    vector<vector<unsigned> >().swap(Labs),

    vector<vector< vector<long long> > >().swap(SearchPath); 

    vector<vector<trielem> >().swap(conR);
    
    vector<unsigned>().swap(w2d); 
    vector<long long>().swap(PathCnt); 

    vector<vector<subt> >().swap(TaskList);
    vector<subt>().swap(subtasks);

    d2w.clear();

    CurW = 0, dmax = 0, tag=0, fflg = 0;
    paths = 0, uness = 0, UpperBound = 0, Wavg = 0;
}


void Parameter(){
    MAXDIS = 2; MAXMOV = 1;

    while( MAXINT / (totalV * 2) >= MAXDIS ) {
        MAXDIS *= 2;
        ++MAXMOV;
    }
    MASK = MAXDIS - 1;
}


void DstValueUpdate(int dis){
    unsigned cur = 0;
    for (int i=0; i<con[dst].size(); ++i){
        unsigned w = con[dst][i] >> MAXMOV, 
                 weight = con[dst][i] & MASK;
        
        if (label[w].size() == 0) continue; // check labels of neighbors

        // get the final element
        unsigned lab = label[w][label[w].size()-1] >> MAXMOV,
                  dd = label[w][label[w].size()-1] & MASK;
        
        if (dd == dis-1){
            unsigned vtrans = lab<weight? lab:weight;
            if (vtrans > cur)
                cur = vtrans; 
        }
    }


    if (CurW < cur){ 
        CurW = cur; // new maximum value
        unsigned elem = CurW << MAXMOV | dis; // Update the label of dst
        label[dst].push_back(elem);
        labelR[src].push_back(elem);
    }

}


void GraphInitial(string filename, int ccc){
    string s;
    const char *filepath = filename.c_str();
    ifstream infile;

    infile.open(filepath);
    if(!infile.is_open()){
        cout<<"No such file!"<<endl;
        exit(-1);
    }

    long xx = 0;

    while(getline(infile, s)){
        char* strc = new char[strlen(s.c_str())+1];
        strcpy(strc, s.c_str());
        char* s1 = strtok(strc," ");
        
        if (xx == 0){
            totalV = atoi(s1);
            con.resize(totalV);
            Parameter();
        }else{
            while(s1){
                int va = xx - 1, vb = atoi(s1)-1,
                    w = (va+vb) % weig + 1;
                // cout<<va<<"  "<<vb<<"  "<<w<<endl;
                
                if (w <= ccc){
                    unsigned val = vb << MAXMOV | w;

                    con[va].push_back(val);
                }
                s1=strtok(NULL," ");
            }
        }

        xx += 1;
        
        delete s1, strc;
    }

    infile.close();

    for( int i = 0; i < totalV; ++i )
		v2degree.push_back(make_pair(con[i].size(), i));

	sort(v2degree.rbegin(), v2degree.rend());
    
}



void GraphInitialWeight(string filename){
    string s;
    const char *filepath = filename.c_str();
    ifstream infile;

    infile.open(filepath);
    if(!infile.is_open()){
        cout<<"No such file!"<<endl;
        exit(-1);
    }

    long xx = 0;

    while(getline(infile, s)){
        char* strc = new char[strlen(s.c_str())+1];
        strcpy(strc, s.c_str());
        char* s1 = strtok(strc," ");
        
        if (xx == 0){
            totalV = atoi(s1);
            con.resize(totalV);
            Parameter();
        }else{
            int va = xx - 1, vb, w;
            int ff = 0;
            while(s1){
                if (ff % 2 == 0)
                    vb = atoi(s1)-1;
                else{
                    w = atoi(s1);
                    // cout<<va<<"  "<<vb<<"  "<<w<<endl;
                    unsigned val = vb << MAXMOV | w;
                    con[va].push_back(val);
                }

                ff += 1;          
                s1=strtok(NULL," ");
            }
        }

        xx += 1;
        
        delete s1, strc;
    }

    infile.close();

    for( int i = 0; i < totalV; ++i )
		v2degree.push_back(make_pair(con[i].size(), i));

	sort(v2degree.rbegin(), v2degree.rend());
    cout<<"size: "<<v2degree.size()<<endl;
}


void LabelMerge(vector<unsigned>& L1, vector<unsigned>& L2, unsigned vid){
    
    for (int i=0; i<L1.size(); ++i){
        
        unsigned w1 = L1[i] >> MAXMOV, d1 = L1[i] & MASK;

        for (int j=0; j<L2.size(); ++j){

            unsigned w2 = L2[j] >> MAXMOV, d2 = L2[j] & MASK;

            if (d1 + d2 > dmax) continue;

            unsigned ww = w1<w2? w1:w2;
            unsigned dd = d1+d2;

            if (d2w[dd] == ww){
                Labs[vid].push_back( (dd<<MAXMOV) | d1 ); // 
            } 
        }
    } 

    if (Labs[vid].size() > 0){
        SearchPath[vid].resize(dmax+1); // record the number of paths in each hop
        for (int ii=0; ii<dmax+1; ++ii){
            SearchPath[vid][ii].resize(dmax+1);
        }
    }
}


unsigned AddEdge(vector<unsigned>& Ls, vector<unsigned>& Ld, unsigned weight){
    
    unsigned val = 0;

    for (int i=0; i<Ls.size(); ++i){
        
        unsigned d1 = Ls[i] >> MAXMOV, ds1 = Ls[i] & MASK, ww = d2w[d1];

        if (weight < ww) break; // 增序，所以后面都是 小于

        for (int j=0; j<Ld.size(); ++j){
            
            unsigned d2 = Ld[j] >> MAXMOV, ds2 = Ld[j] & MASK;
            
            if (d1 != d2) continue; // 在同一个距离的路径上

            if (ds1 + 1 == ds2){
                val += (1 << d1);
            }
        }
    }

    return val;
}


void ParallelBuild(){
    
    omp_set_num_threads(threads);
    
    label.resize(totalV), labelR.resize(totalV);
    
    vector<unsigned> Tv1(totalV, 0), Tv2(totalV, 0);

    label[src].push_back(((weig+1) << MAXMOV | 0)); // weight+dis
    
    labelR[dst].push_back(((weig+1) << MAXMOV | 0));

    int dis = 1;

    for( long long cnt1 = 1, cnt2 = 1; dis <= MAXDIS; ++dis ){
    
    cnt1 = 0, cnt2 = 0;
    
    DstValueUpdate(dis); // Confirm CurW of the dis-path

    # pragma omp parallel
    {
        int pid = omp_get_thread_num(), np = omp_get_num_threads();

        for( int u = pid; u < totalV; u += np ){
            
            if (u == src or u == dst) continue; // src 不需要考虑，dst是优先计算

            unsigned Th1 = 0, Th2 = 0, Beta1 = CurW, Beta2 = CurW;
            
            if (label[u].size() > 0){
                unsigned val = label[u][label[u].size()-1] >> MAXMOV;
                if (val > Beta1) Beta1 = val;
            }

            if (labelR[u].size() > 0){
                unsigned val = labelR[u][labelR[u].size()-1] >> MAXMOV;
                if (val > Beta2) 
                    Beta2 = val;
            }

            // ====================================
            for( int i = 0; i < con[u].size(); ++i ){
            
                unsigned w = con[u][i] >> MAXMOV, weight = con[u][i] & MASK;
                
                if ((w == src or w == dst) and dis >= 2) 
                    continue; // 源于 src / dst 的信息不用考虑

                if (label[w].size() > 0){
                    // get the final element
                    unsigned lab = label[w][label[w].size()-1] >> MAXMOV,
                            dd = label[w][label[w].size()-1] & MASK;
                    
                    if (dd == dis-1){
                        unsigned vtrans = lab<weight? lab:weight;
                        Th1 = vtrans>Th1? vtrans:Th1;
                    }
                }

                if (labelR[w].size() > 0){
                    unsigned lab = labelR[w][labelR[w].size()-1] >> MAXMOV,
                            dd = labelR[w][labelR[w].size()-1] & MASK;

                    if (dd == dis-1){
                        unsigned vtrans = lab<weight? lab:weight;
                        Th2 = vtrans>Th2? vtrans:Th2;
                    }
                }
            }

            // == Thresold v.s. beta ==
            if (Th1 > Beta1){
                #pragma omp critical
                {
                    cnt1 += 1;
                }
                Tv1[u] = Th1 << MAXMOV | dis;
            }  

            if (Th2 > Beta2){
                #pragma omp critical
                {
                    cnt2 += 1;
                }
                Tv2[u] = Th2 << MAXMOV | dis;
            }  
        }

    }

    # pragma omp parallel
    {
        int pid = omp_get_thread_num(), np = omp_get_num_threads();
        
        for( int u = pid; u < totalV; u += np ){
            if (Tv1[u] != 0){
                label[u].push_back(Tv1[u]);
                Tv1[u] = 0;
            }

            if (Tv2[u] != 0){
                labelR[u].push_back(Tv2[u]);
                Tv2[u] = 0;
            }
        }
    }

    
    // cout<<"dis: "<<dis<<"  Thre: "<<CurW<<"  "<<cnt1<<"  "<<cnt2<<endl;

    if (cnt1 == 0 or cnt2 == 0) break;
    }

}


void GraphBuild(){
    // == Acticated from dst ==

    Labs.resize(totalV);
    SearchPath.resize(totalV);

    unsigned lst = label[dst].size()-1;
    
    dmax = label[dst][lst] & MASK;
    w2d.resize(weig+1);

    if (dmax > 0){
    for (int i=0; i<label[dst].size(); ++i){
        unsigned lab = label[dst][i] >> MAXMOV, dd = label[dst][i] & MASK;
        d2w[dd] = lab;
        tag += (1<<dd);

        Labs[src].push_back( ((dd<<MAXMOV)|0) );  
        Labs[dst].push_back( ((dd<<MAXMOV)|dd) );
        
        if (i > 0){
            unsigned st_lab = label[dst][i-1] >> MAXMOV,
                    st_dis = label[dst][i-1] & MASK;
                    
            for (int j=st_lab; j<lab; ++j)
                w2d[j] = st_dis;
        }

        if ( i == label[dst].size()-1){
            for (int iii=lab; iii<=weig; ++iii)
                w2d[iii] = dd;
        }
    }

    omp_set_num_threads(threads);
    # pragma omp parallel
    {
        int pid = omp_get_thread_num(), np = omp_get_num_threads();

        for( int u = pid; u < totalV; u += np ){
            
            if (u == src or u == dst) continue;

            vector<unsigned>& L1 = label[u];
            vector<unsigned>& L2 = labelR[u];

            if (L1.size() == 0 or L2.size() == 0) continue;

            LabelMerge(L1, L2, u);
        }
    }

    conR.resize(totalV);
    omp_set_num_threads(threads);
    # pragma omp parallel
    {
        int pid = omp_get_thread_num(), np = omp_get_num_threads();

        for( int u = pid; u < totalV; u += np ){

            if (Labs[u].size() == 0) continue;

            vector<unsigned>& adj = con[u];

            for (int i=0; i<adj.size(); ++i){

                unsigned vid = adj[i] >> MAXMOV, weight = adj[i] & MASK;

                if (Labs[vid].size() == 0) continue;

                if (u == vid) continue;

                unsigned val = AddEdge(Labs[u], Labs[vid], weight);

                if (val > 0){
                    trielem e1(vid, weight, val);
                    conR[u].push_back(e1); // 不需要权重了
                }
            }
        
            sort(conR[u].begin(), conR[u].end(), cmp);
        }
    }

    }
    
}


void GraphBuild_Single(){
    // == Acticated from dst ==

    Labs.resize(totalV);
    SearchPath.resize(totalV);

    int lst = label[dst].size()-1;

    if (lst >= 0){
        dmax = label[dst][lst] & MASK;
        w2d.resize(weig+1);
        for (int i=0; i<label[dst].size(); ++i){
            unsigned lab = label[dst][i] >> MAXMOV, dd = label[dst][i] & MASK;
            d2w[dd] = lab;
            tag += (1<<dd);

            Labs[src].push_back( ((dd<<MAXMOV)|0) );  
            Labs[dst].push_back( ((dd<<MAXMOV)|dd) );
            
            if (i > 0){
                unsigned st_lab = label[dst][i-1] >> MAXMOV,
                        st_dis = label[dst][i-1] & MASK;
                        
                for (int j=st_lab; j<lab; ++j)
                    w2d[j] = st_dis;
            }

            if ( i == label[dst].size()-1){
                for (int iii=lab; iii<=weig; ++iii)
                    w2d[iii] = dd;
            }
        }


        for( int vid = 0; vid < totalV; vid ++ ){
            
            if (vid == src or vid == dst) continue;

            vector<unsigned>& L1 = label[vid];
            vector<unsigned>& L2 = labelR[vid];

            if (L1.size() == 0 or L2.size() == 0) continue;

            LabelMerge(L1, L2, vid);
        }
        
        conR.resize(totalV);

        for( int vv = 0; vv < totalV; vv ++ ){

            if (Labs[vv].size() == 0) continue;

            vector<unsigned>& adj = con[vv];

            for (int i=0; i<adj.size(); ++i){

                unsigned vid = adj[i] >> MAXMOV, weight = adj[i] & MASK;

                if (Labs[vid].size() == 0) continue;

                if (vv == vid) continue;

                unsigned val = AddEdge(Labs[vv], Labs[vid], weight);

                if (val > 0){
                    trielem e1(vid, weight, val);
                    conR[vv].push_back(e1); // 不需要权重了
                }
            }
        
            sort(conR[vv].begin(), conR[vv].end(), cmp);
        }
    }

    


    
}


void DFSearch(unsigned u, unsigned tgt,
              vector<unsigned>& s_, vector<unsigned>& lab_){
    
    s_.push_back(u);
    unsigned tag1  = lab_[lab_.size()-1] >> MAXMOV, 
             wmin = lab_[lab_.size()-1] & MASK;

    for (int i=0; i<conR[u].size(); ++i){

        unsigned v = conR[u][i].vid, w = conR[u][i].wval, 
                 t = conR[u][i].flgv;

        unsigned wmm = wmin < w? wmin:w,
                 dd = w2d[wmm];

        unsigned newTag = tag1 & t;

        if (newTag == 0) continue;

        if (v == dst){
            unsigned valu = (1<<(s_.size())), flg = valu & newTag;
            if (flg > 0)  
                paths += 1;

        }else if(s_.size() == dd-1){
            paths += 1;

        }else if (s_.size() < dd){
            unsigned vval = newTag<<MAXMOV | wmm;
            lab_.push_back(vval);

            DFSearch(v, tgt, s_, lab_);

            lab_.pop_back();

        }else{
            uness += 1;
        }
    }
    
	s_.pop_back();

    return;
}


void DFSearchParallel(unsigned u, unsigned tgt,
              vector<unsigned>& s_, vector<unsigned>& lab_,
              long long& pcnt){
    
    s_.push_back(u);
    unsigned tag1 = lab_[lab_.size()-1] >> MAXMOV, 
             wmin = lab_[lab_.size()-1] & MASK;

    for (int i=0; i<conR[u].size(); ++i){

        unsigned v = conR[u][i].vid, w = conR[u][i].wval, 
                 t = conR[u][i].flgv;

        unsigned wmm = wmin < w? wmin:w,
                 dd = w2d[wmm];

        unsigned newTag = tag1 & t;

        if (newTag == 0) continue;

        if (v == dst){
            unsigned valu = (1<<(s_.size())), flg = valu & newTag;
            if (flg > 0) pcnt += 1;
        }else if(s_.size() == dd-1){
            pcnt += 1;
        }else if(s_.size() > dd){
            break; // 后续的值肯定对应dis更小
        }else if (s_.size() < dd){
            unsigned vval = newTag<<MAXMOV | wmm;
            lab_.push_back(vval);

            DFSearchParallel(v, tgt, s_, lab_, pcnt);

            lab_.pop_back();

        }else{
            uness += 1;
        }
    }
    
	s_.pop_back();

    return;
}


void PathUpdate(unsigned v1, unsigned v2, unsigned w, unsigned k){ // v2 is a neighbor of v1
    vector<unsigned>& L1 = Labs[v1];
    vector<unsigned>& L2 = Labs[v2];

    for (int i=0; i<L1.size(); ++i){
        
        unsigned d1 = L1[i] >> MAXMOV, ds1 = L1[i] & MASK, wmin = d2w[d1];

        if (w < wmin) break;

        for (int j=0; j<L2.size(); ++j){
            
            unsigned d2 = L2[j] >> MAXMOV, ds2 = L2[j] & MASK;

            if (d1 == d2 and ds2-ds1 == 1){
                // cout<<"V: "<<v1<<"  "<<v2<<"  hop: "<<k<<"   distance: "<<d1<<endl;
                SearchPath[v1][d1][k] += SearchPath[v2][d1][k-1];
            }
        }
    }
}


void PredictR(){

    SearchPath[src].resize(dmax+1);
    SearchPath[dst].resize(dmax+1);

    for (int i=0; i<dmax+1; ++i){
        SearchPath[src][i].resize(dmax+1);
        SearchPath[dst][i].resize(dmax+1);
    }

    for (int i=0; i<Labs[dst].size(); ++i){
        unsigned d = Labs[dst][i] >> MAXMOV, ds = Labs[dst][i] & MASK;
        SearchPath[dst][d][d-ds] = 1; 
    }

    for (int i=1; i<dmax+1; ++i){

        for (int j=0; j<conR.size(); ++j){

            vector<trielem>& neigs = conR[j];

            if (neigs.size() == 0) 
                continue;

            for (int k=0; k<neigs.size(); ++k){

                unsigned vid = neigs[k].vid, w = neigs[k].wval;

                PathUpdate(j, vid, w, i);
            }
        }
    }

    
    for (int i=0; i<dmax+1; ++i){

        if (d2w[i] == 0) continue;

        UpperBound += SearchPath[src][i][i];
        // cout<<"d: "<<i<<"   Paths: "<<SearchPath[src][i][i]<<endl;
        if (fflg == 0 and UpperBound > MaxLoad){
            fflg = 1;
        }
    }

    Wavg = UpperBound / threads + 1;
}


void TaskDivision(unsigned u, unsigned tgt,
              vector<unsigned>& s_, 
              vector<unsigned>& lab_){
    
    s_.push_back(u);           
    unsigned tag1  = lab_[lab_.size()-1] >> MAXMOV, 
             wmin = lab_[lab_.size()-1] & MASK;

    for (int i=0; i<conR[u].size(); ++i){
        unsigned v = conR[u][i].vid, w = conR[u][i].wval, 
                 t = conR[u][i].flgv;

        unsigned wmm = wmin < w? wmin:w,
                 dd = w2d[wmm];

        unsigned newTag = tag1 & t;

        if (newTag == 0) continue;

        if (v == dst){
            unsigned valu = (1<<(s_.size())), flg = valu & newTag;
            if (flg > 0) paths += 1;

        }else if(s_.size() == dd-1){
            paths += 1;

        }else if (s_.size() < dd){
            vector< vector<long long> >& slist = SearchPath[v];
            long long workload = 0;
            int d = s_.size();

            for (int kkk=0; kkk<Labs[v].size(); ++kkk){
                unsigned dval = Labs[v][kkk] >> MAXMOV;

                if (dval-d > 0 and dval <= dd) {
                    workload += slist[dval][dval-d];
                }
            }

            // cout<<u<<" - "<<v<<" * "<<Labs[v].size()<<"   "<<workload<<endl;

            unsigned vval = newTag<<MAXMOV | wmm;
            
            lab_.push_back(vval);

            if (workload <= Wavg){
                s_.push_back(v);
                subt elem(workload, s_, lab_);              
                subtasks.push_back(elem);
                s_.pop_back();
            }else
                TaskDivision(v, tgt, s_, lab_);

            lab_.pop_back();
        }else{
            uness += 1;
        }

    }

    s_.pop_back();
}


void QueryPairs(string ff){
    string filename = "Result/small/"+ff+"_result";
    string s;
    const char *filepath = filename.c_str();
    ifstream infile;

    infile.open(filepath);
    if(!infile.is_open()){
        cout<<"No such file!"<<endl;
        exit(-1);
    }

    long paths;
    int vs, vt, dmax;
    float t1, t2;

    while(getline(infile, s)){
        char* strc = new char[strlen(s.c_str())+1];
        strcpy(strc, s.c_str());
        // cout<<strc<<endl;
        char* s1 = strtok(strc," ");
        vs = atoi(s1);
        
        int flg = 0;

        while(s1){
            
            s1=strtok(NULL," ");
            
            switch (flg){
                case 0: vt = atoi(s1); break;
                case 1: dmax = atoi(s1); break;
                case 2: t1 = atof(s1); break;
                case 3: t2 = atof(s1); break;
                case 4: paths = atol(s1); break;
            }

            if (flg == 4)
                QPairs.push_back(make_pair(vs,vt));

            flg += 1;
        }
        
        delete s1, strc;
    }

    infile.close();
}

void VertexDivide(){
    vector<long long> Paths(threads);

    for (int i=0; i<subtasks.size(); ++i){
        int minP = min_element(Paths.begin(),Paths.end()) - Paths.begin();
        Paths[minP] += subtasks[i].cnts;
        TaskList[minP].push_back(subtasks[i]);
    }

    for (int i=0; i<Paths.size(); ++i){
        cout<<Paths[i]<<endl;
    }

}

int main(){

    int ccc = 1000000;
    weig = 10;

    string na = "ukunion";
    string name = na+".graph";
    string Newfile = na+"_result", Newfile2 = na+"_Query";

    threads = 10;
    vector<unsigned> stk, lab;
    MaxLoad = 100000000;
    int TaskCnt = 10000;

    srand((unsigned)time(NULL));
    
    GraphInitial(name, ccc);

	const char *file = Newfile.c_str(), *file1 = Newfile2.c_str();
	fstream outfileX, outfileX1 ;
	outfileX.open(file, ios::app);
	outfileX1.open(file1, ios::app);

    for (int ij=0; ij <TaskCnt; ++ij){

        aplace = rand() % 1000000, bplace = rand() % (totalV-1000000) + 1000000;
        src = v2degree[aplace].second,  dst = v2degree[bplace].second;

        float t = omp_get_wtime();

        ParallelBuild();

        float t1 = omp_get_wtime();

        GraphBuild_Single();

        if (dmax == 0){
            InitialP();
            continue;
        }
 
        PredictR();

        if (fflg == 1){
            InitialP();
            continue;
        } 

        if (UpperBound > MaxLoad){
            string inf = to_string(src) + " " + to_string(dst) + " " + to_string(dmax)+ " " +to_string(UpperBound);
            outfileX1<<inf<<endl;
            InitialP();
            continue;
        }

        float t2 = omp_get_wtime();
        
        lab.push_back(tag<<MAXMOV | (weig+1));

        if (UpperBound < MaxLoad){ // slight workload
            DFSearch(src, dst, stk, lab);
        }else{
            TaskDivision(src, dst, stk, lab);

            TaskList.resize(threads);
            
            sort(subtasks.begin(), subtasks.end(), cmp1);

            VertexDivide();

            # pragma omp parallel
            {
                int pid = omp_get_thread_num(), np = omp_get_num_threads();
                
                for( int u = pid; u < TaskList.size(); u += np ){
                    
                    long long cntt = 0;
                    
                    for (int kk=0; kk<TaskList[u].size(); ++kk){
                        
                        vector<unsigned>& tsk = TaskList[u][kk].que;
                        vector<unsigned>& labb = TaskList[u][kk].lab;
                        
                        unsigned vid = tsk[tsk.size()-1];
                        tsk.pop_back();

                        DFSearchParallel(vid, dst, tsk, labb, cntt);
                    }

                    #pragma omp critical
                    {
                        paths += cntt;
                    }
                }
            }
        }
        lab.clear();

        float t3 = omp_get_wtime();

        cout<<src<<" * "<<dst<<" * "<<t3-t<<" * "<<paths<<endl;

        string inf = to_string(src) + " " + to_string(dst) + " " +
                     to_string(dmax) + " " + to_string(t3-t) + " " + to_string(t3-t2) + " " + to_string(paths); 
        
        outfileX<<inf<<endl;

        InitialP();  
    }


    return 0;
}
