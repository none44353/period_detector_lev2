#ifndef _detectorHG_H
#define _detectorHG_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstring>
#include "BOBHash32.h"
#include "BOBHash64.h"

using namespace std;

#define T 20
#define MAXM 500005
#define MAXL 1000005 

class AlgHG //基于Heavy Guardian的版本
{   
  private:
    int t, M;
    BOBHash32* bobhash;

    double var[T * MAXM];
    uint64_t id[T * MAXM];
    int counter[T * MAXM]; 
  
  public:
    AlgHG(int t, int M) : t(t), M(M) { 
        bobhash = new BOBHash32(123);
        memset(counter, 0, sizeof(counter));
    } 
    
    /*
    bool check_near(double c, double x) {
        return fabs(x - c) <= delta;
    }*/
    bool check_near_n(double c, double x) {
      //  return fabs(x - c) <= delta;
         return fabs(x - c) <= 2 * delta;
    }

    void insert(const uint64_t& s, const double& x) { //丢一个元素进来
        if (x < rangeL) return;
        if (x > rangeR) return;
        unsigned long long int H = bobhash -> run((char *)&s, 8);
        unsigned long long int pos = H % M;
        
        int nl = pos * t, nr = (pos + 1) * t, cnt;
        
        for (int i = nl; i < nr; ++i) if ((cnt = counter[i]) && id[i] == s) {
            double nvar = var[i] * cnt / (cnt + 1) + x / (cnt + 1);
            if (check_near_n(var[i], x)) {
                var[i] = nvar;
                ++counter[i];
            }
            else {
                int p = rand() % (counter[i] + 1);
                if (!p) { //替换
                    var[i] = x;
                   // counter[i] = 1; 
                    ++counter[i];
                }
            }
            return;
        } 
        
        int nx = nl;
        for (int i = nl; i < nr; ++i) 
            if (counter[i] < counter[nx]) nx = i;
        
        int p = rand() % (counter[nx] + 1);
        if (!p) { //替换
            id[nx] = s;
            var[nx] = x;
            counter[nx] = 1; 
        }

        return;
    }

    pair <int, double> query(const uint64_t& s) { //询问某个id是否有周期性
        unsigned long long int H = bobhash -> run((char *)&s, 8);
        unsigned long long int pos = H % M;
        
        int nl = pos * t, nr = (pos + 1) * t, cnt;
        
        for (int i = nl; i < nr; ++i) if ((cnt = counter[i]) && id[i] == s) 
            return make_pair(cnt, var[i]);
        
        return make_pair(0, -1);
    }
};
#endif
