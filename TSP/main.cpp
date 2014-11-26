//
//  main.cpp
//  TSP
//
//  Created by Sami Purmonen on 25/11/14.
//  Copyright (c) 2014 Sami Purmonen. All rights reserved.
//


#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <map>
#include <assert.h>

using namespace std;


class TravelingSalesmanProblem{
public:
    int points;
    double *x;
    double *y;
    
    static TravelingSalesmanProblem *createFromStdin(){
        TravelingSalesmanProblem *instance = new TravelingSalesmanProblem();;
        cin >> instance->points;
        instance->y = new double[instance->points];
        instance->x = new double[instance->points];
        for (auto i = 0; i < instance->points; i++)
            cin >> instance->x[i] >> instance->y[i];
        return instance;
    }
    
    static TravelingSalesmanProblem *testInstance(){
        TravelingSalesmanProblem *instance = new TravelingSalesmanProblem();;
        instance->points = 10;
        instance->x = new double[instance->points];
        instance->y = new double[instance->points];
        instance->x[0] = 95.0129;		//	95.0129 61.5432
        instance->y[0] = 61.5432;
        instance->x[1] = 23.1139;		//	23.1139 79.1937
        instance->y[1] = 79.1937;
        instance->x[2] = 60.6843;		//	60.6843 92.1813
        instance->y[2] = 92.1813;
        instance->x[3] = 48.5982;		//	48.5982 73.8207
        instance->y[3] = 73.8207;
        instance->x[4] = 89.1299;		//	89.1299 17.6266
        instance->y[4] = 17.6266;
        instance->x[5] = 76.2097;		//	76.2097 40.5706
        instance->y[5] = 40.5706;
        instance->x[6] = 45.6468;		//	45.6468 93.5470
        instance->y[6] = 93.5470;
        instance->x[7] = 1.8504;		//	1.8504 91.6904
        instance->y[7] = 91.6904;
        instance->x[8] = 82.1407;		//	82.1407 41.0270
        instance->y[8] = 41.0270;
        instance->x[9] = 44.4703;		//	44.4703 89.3650
        instance->y[9] = 89.3650;
        return instance;
    }
    
    double dist(int a, int b){
        return sqrt((x[a] - x[b])*(x[a] - x[b]) + (y[a] - y[b])*(y[a] - y[b]));
    }
    
    double best = 9999999;
    void solvenaiveExactly(int pointsLeft, double distance, int current){
        if (pointsLeft == 1) best = min(best, distance + dist(current, 0));
        
        if (distance >= best) return;
        for (auto i = 1; i < pointsLeft; i++)
        {
            swap(x[i], x[pointsLeft - 1]);
            swap(y[i], y[pointsLeft - 1]);
            solvenaiveExactly(pointsLeft - 1, distance + dist(current, pointsLeft - 1), pointsLeft - 1);
            swap(y[i], y[pointsLeft - 1]);
            swap(x[i], x[pointsLeft - 1]);
        }
    }
    
    void solvenaive2(int pointsLeft, double distance, int current){
        if (pointsLeft == 1) best = min(best, distance + dist(current, 0));
        
        
        auto bestIndex = 0;
        double closestDistance = 9999999;
        for (auto i = 1; i < pointsLeft; i++)
        {
            if (dist(current, i) < closestDistance){
                closestDistance = dist(current, i);
                bestIndex = i;
            }
        }
        swap(x[bestIndex], x[pointsLeft - 1]);
        swap(y[bestIndex], y[pointsLeft - 1]);
        solvenaiveExactly(pointsLeft - 1, distance + closestDistance, pointsLeft - 1);
        swap(y[bestIndex], y[pointsLeft - 1]);
        swap(x[bestIndex], x[pointsLeft - 1]);
        
    }
    
    void solvenaive(vector<int> &used, double distance, int current, int depth){
        if (used.size() == points) {
            best = distance + dist(current, 0);
            return;
        }
        
        auto bestIndex = 0;
        double closestDistance = 9999999;
        for (int i = 0; i < points; i++){
            if (find(used.begin(), used.end(), i) != used.end()) continue;
            if (dist(current, i) < closestDistance){
                closestDistance = dist(current, i);
                bestIndex = i;
            }
        }
        
        used.push_back(bestIndex);
        solvenaive(used, distance + closestDistance, bestIndex, depth + 1);
    }
    
    void printBits(unsigned long bits) {
        for (int i = 31; i >= 0; i--) {
            cout << ((bits >> i) & 1) << " ";
        }
        cout << endl;
    }
    
    double dynamicExact() {
        int n = points;
        map<unsigned long, map<int, int>> parent;
        map<unsigned long, map<int, double>> C;
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                unsigned long v = 0;
                v |= 1 << j;
                v |= 1 << k;
                //printBits(v);
                C[v][k] = dist(j, k);
            }
        }
        for (int s = 3; s <= n; s++) {
            // Stack 1,2,3
            unsigned long S = 0;
            for (int i = 0; i < s; i++) {
                S |= 1 << i;
            }
            //            printBits(S);
            
            /* Bit rage */
            int currentBit = 0;
            while (true) {
                //cout << "s = " << S << endl;
                for (int k = 1; k < n; k++) {
                    if ((S >> k) & 1) {
                        double minimum = 1e300;
                        for (int m = 1; m < n; m++) {
                            if (m != k && ((1<<m) & S)) {
                                // k = nuvarande slutplats
                                // m = platsen innan nuvarande slutplats
                                unsigned long tmpS = S ^ (1 << k);
                                double distance = dist(m, k);
                                //cout << "S = " << S << endl;
                                //cout << "tmpS = " << tmpS << endl;
                                //cout << "m = " << m << endl;
                                //cout << "k = " << k << endl;
//                                assert(C[tmpS][m] != 0);
                                //8 4 1
                                //1101 1
                                //1111
                                double tmpVal = C[tmpS][m] + distance;
                                if (tmpVal < minimum) {
                                    parent[S][k] = m;
                                    C[S][k] = tmpVal;
                                    minimum = tmpVal;
                                }
                            }
                        }
                    }
                }
                if (currentBit >= n - s) {
                    break;
                }
                if (((S >> (currentBit + 1)) & 1) == 0) {
                    S += (1 << currentBit);
                    currentBit++;
                }
                else {
                    int collidingOnes = 1;
                    int i = 2;
                    while ((S >> (currentBit + i)) & 1) {
                        collidingOnes++;
                        i++;
                    }
                    S += (1 << currentBit);
                    for (int i = 0; i < collidingOnes; i++) {
                        S |= 1 << i;
                    }
                    currentBit = 0;
                }
            }
        }
        
        double opt = 9999999;
        //        printBits((1 << ()n) - 1);
        int node = -1;
        unsigned long S = (1 << n) - 1;
        for (int k = 1; k < n; k++) {
            double tmp = C[S][k];
            opt = min(tmp, opt);
            node = k;
        }

        while (true) {
            cout << node << endl;
            if (node == 0) {
                break;
            }
            int oldNode = node;
            node = parent[S][node];
            S = S^(1 << oldNode);
        }
        return opt;
    }
    
    vector<int> greedy() {
        vector<int> vector;
        int *tour = new int[points];
        bool *used = new bool[points];
        for (int i = 0; i < points; i++) {
            used[i] = false;
        }
        tour[0] = 0;
        used[0] = true;
        for (int i = 1; i < points; i++) {
            int best = -1;
            for (int j = 0; j < points; j++) {
                if (!used[j] && (best == -1 || dist(tour[i - 1], j) < dist(tour[i - 1], best))) {
                    best = j;
                }
            }
            tour[i] = best;
            used[best] = true;
        }
        for (int i = 0; i < points; i++) {
            vector.push_back(tour[i]);
        }
        return vector;
    }
};


int main(int argc, char **argv){
    auto instance = TravelingSalesmanProblem::createFromStdin();
    instance->dynamicExact();
    
    //    for (auto node: instance->greedy()) {
    //        cout << node << endl;
    //    }
    //    cout << instance->best << endl;
    
    //    instance->solvenaiveExactly(10, 0, 0);
    //    cout << instance->best << endl;
    
    
    delete instance;
    return 0;
}