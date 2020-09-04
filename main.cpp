#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <time.h>

using namespace std;

#define Range 500
#define delta 0.01
#define rangeL 0.05
#define rangeR 0.30


#include "ssummary.h"
#include "BF.h"
//#include "detector.h"
#include "detectorHG.h"

//string datapath[60] = {"./130000.dat"};


string datapath[60] = {"../../usr/share/dataset/CAIDA2018/dataset/130000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/130900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/131900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/132900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/133900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134600.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/134900.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135000.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135100.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135200.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135300.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135400.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135500.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135600.dat", 
                "../../usr/share/dataset/CAIDA2018/dataset/135700.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135800.dat",
                "../../usr/share/dataset/CAIDA2018/dataset/135900.dat"};

ifstream fin;

const int M = 2e7;
const int _K = 500;

pair <uint64_t, double> Read()//新CAIDA
{   
    static bool isfirstRead = true;
    static int curFinID = 0;
    static double offset = 0;
    static double lastT = 0;

	double t; uint64_t s, _s;
    if (isfirstRead) {
        isfirstRead = false;
        fin.open(datapath[curFinID], std :: ios :: binary);
    }

	if (fin.eof()) {
        fin.close();
        fin.open(datapath[++curFinID], std :: ios :: binary);
        if (curFinID > 60) {
            fin.close();
            exit(0);
        }
    }
    fin.read((char*)&s, sizeof(uint64_t)); //srcip(4)+dstip(4) 
    fin.read((char*)&_s, 5);//srcport destport protcol

    fin.read((char*)&t, sizeof(double));

    
    t += offset;
    if(t < lastT) {
        offset += lastT - t;
        t += (lastT - t);
    }
    lastT = t;

	return make_pair(s, t * 2); 
}

pair <uint64_t, double> input[M + 7];

bool check_near(double c, double x) {
    return fabs(x - c) <= delta;
}

struct getGT{
    map <uint64_t, int> count;
    map <uint64_t, set <double> > table;
    
    void init() {
        count.clear();
        table.clear();
    }

    void insert(uint64_t id, double key) {
        if (key < rangeL) return;
        if (key > rangeR) return;
        if (count.find(id) == count.end()) 
            count[id] = 0, table[id].clear();
        count[id]++, table[id].insert(key);
    }

    int query(const uint64_t& id) {
        if (count.find(id) == count.end()) return 0;
        int n = count[id], mx = 0, cur = 1; 
        
        set <double> :: iterator itl = table[id].begin();
        set <double> :: iterator itr = table[id].begin();

        for (; itr != table[id].end(); ++itr, ++cur) {
            while (*itl < *itr - 2 * delta) ++itl, --cur;
            mx = max(mx, cur);
        }

        return mx;
    }

    int querynumber(const uint64_t& id, double c) {
        set <double> :: iterator it = table[id].begin();
        int cnt = 0;
        for (; it != table[id].end(); ++it)
            if (check_near(c, *it)) ++cnt;
        return cnt;   
    }

}intervalGT;


map <uint64_t, double> timeStamp;
map <uint64_t, int> intervalAnswer;
map <uint64_t, bool> check;

int total_ele, total_save;
pair <int, uint64_t> ele[M + 7], save[M + 7];
set <uint64_t> topKid, topKid_ours;


void GroundTruth() {
    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 

        uint64_t id = e.first;
        //if (id != 17019272071912997107ULL) continue;
        //if (id != 17213472216328176400ULL) continue;
        // if (id != 6087492077255636877ULL) continue;
        double curTime = e.second;
        double lastTime = -1;
        
        if (timeStamp.find(id) != timeStamp.end()) {
            lastTime = timeStamp[id];
         //  if (curTime - lastTime > rangeL) printf("#%.6lf\n", curTime - lastTime);
            intervalGT.insert(id, curTime - lastTime);
        }
        
        timeStamp[id] = curTime;
    }

    intervalAnswer.clear();
    total_ele = 0;
    for (int i = 0; i < M; ++i){
        auto e = input[i]; 
        uint64_t id = e.first;
        
       // if (id != 17019272071912997107ULL) continue;
      ////  if (id != 17213472216328176400ULL) continue;
      //  if (id != 6087492077255636877ULL) continue;
        if (intervalAnswer.find(id) != intervalAnswer.end()) continue;
        //printf("")
        int result = intervalGT.query(id);
        intervalAnswer[id] = result;
        ele[++total_ele] = make_pair(result, id);
      //  printf("c %.6lf per %.6lf\n", result.second, intervalGT.queryPercentage(id, result.second));
    }

    sort(ele + 1, ele + total_ele + 1);
    for (int k = total_ele, i = 1; i <= _K; ++i, --k)
        topKid.insert(ele[k].second); //printf("top%d %llu %d\n", total_ele - k + 1, (long long unsigned int)ele[k].second, intervalAnswer[ele[k].second]);
}


int precision[2], recall[2];
int cnt;

int meanprecision[2];

void OURS() {
    bloomfliter* ntimeStamp = new bloomfliter(4, 32768); //t * M * 8 / 1024 = 1024KB
    AlgHG* HGDetector = new AlgHG(4, 800);  //L * 24 / 1024 = 150KB

    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 

        uint64_t id = e.first;
        double curTime = e.second;
        double lastTime = ntimeStamp -> query(id);
        
        HGDetector -> insert(id, curTime - lastTime);

        ntimeStamp -> insert(id, curTime);
    }


//if (id != 2051345630196722313ULL) continue;
// if (id != 17213472216328176400ULL) continue;
//  if (id != 17019272071912997107ULL) continue;
//#409008156058700588
//#3503478137803909614
//#6087492077255636877
    check.clear();
    total_save = 0;
    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 
        uint64_t id = e.first;
        
        if (check.find(id) != check.end()) continue;
        check[id] = true; ++cnt;
        pair <int, double> guess = HGDetector -> query(id);
        if (guess.first) save[++total_save] = make_pair(guess.first, id);//topKid_ours.insert(id);

       // if (topKid.find(id) != topKid.end()) {
        if (guess.first) {
            ++meanprecision[1];
            int num = intervalGT.querynumber(id, guess.second);
            if (num >= 0.9 * guess.first) ++meanprecision[0];
        }
    }

    sort(save + 1, save + total_save + 1);
    for (int k = total_save, i = 1; i <= total_save; ++i, --k)
        topKid_ours.insert(save[k].second);

    for (auto it = topKid.begin(); it != topKid.end(); ++it) {
        ++recall[1];
        if (topKid_ours.find(*it) != topKid_ours.end()) ++recall[0];
    }

    for (auto it = topKid_ours.begin(); it != topKid_ours.end(); ++it) {
        ++precision[1];
        if (topKid.find(*it) != topKid.end()) ++precision[0];
    }
}

//17213472216328176400 
//1006749987763913488 
int main() {
    srand(0);
    for (int i = 0; i < M + 1; ++i) input[i] = Read();
    
    GroundTruth(); puts("calc GT");

    OURS(); puts("calc OURS");
    printf(" %d %d %.6lf\n", meanprecision[1], meanprecision[0], (double)meanprecision[0]/meanprecision[1]);
    printf("%d\n", cnt);
    printf("#%d %d\n", total_ele, total_save);
    printf("precision %d %d %.6lf\n", precision[1], precision[0], (double)precision[0] / precision[1]);
    printf("recall %d %d %.6lf\n", recall[1], recall[0], (double)recall[0] / recall[1]);
   // printf("centre_precision %.lf %.6lf\n", (double)centre_precision[1], (double)(centre_precision[0] / centre_precision[1]));

	return 0;
}