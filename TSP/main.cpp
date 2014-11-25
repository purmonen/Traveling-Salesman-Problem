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

using namespace std;



class TravelingSalesmanProblem{
public:
    int points;
    
    double *x;
    double *y;
    
    static TravelingSalesmanProblem *createFromStdin(){
        TravelingSalesmanProblem *instance = new TravelingSalesmanProblem();;
        cin >> instance->points;
        instance->x = new double[instance->points];
        instance->y = new double[instance->points];
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
    void solveNaieveExactly(int pointsLeft, double distance, int current){
        if (pointsLeft == 1) best = min(best, distance + dist(current, 0));
        
        if (distance >= best) return;
        for (auto i = 1; i < pointsLeft; i++)
        {
            swap(x[i], x[pointsLeft - 1]);
            swap(y[i], y[pointsLeft - 1]);
            solveNaieveExactly(pointsLeft - 1, distance + dist(current, pointsLeft - 1), pointsLeft - 1);
            swap(y[i], y[pointsLeft - 1]);
            swap(x[i], x[pointsLeft - 1]);
        }
        
    }
    
    void solveNaieve2(int pointsLeft, double distance, int current){
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
        solveNaieveExactly(pointsLeft - 1, distance + closestDistance, pointsLeft - 1);
        swap(y[bestIndex], y[pointsLeft - 1]);
        swap(x[bestIndex], x[pointsLeft - 1]);
        
    }
    
    void solveNaieve(vector<int> used, double distance, int current, int depth){
        if (used.size() == points) { best = distance + dist(current, 0); return; }
        
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
        solveNaieve(used, distance + closestDistance, bestIndex, depth + 1);
        
    }
    
};


int main(int argc, char **argv){
    
    
    auto instance = TravelingSalesmanProblem::testInstance();
    vector<int> vector;
    instance->solveNaieve(vector, 0, 0, 0);
    vector.clear();
    cout << instance->best << endl;
    
    instance->solveNaieveExactly(10, 0, 0);
    cout << instance->best << endl;
    delete instance;
    return 0;
}