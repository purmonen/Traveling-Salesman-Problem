#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <chrono>
#include <fstream>
#include <map>

#ifndef DEBUG
#define assert void
#endif

using namespace std;
auto startTime = chrono::high_resolution_clock::now();

inline bool timeIsRunningOut() {
#ifdef DEBUG
    return startTime + chrono::milliseconds(100000) < chrono::high_resolution_clock::now();
#else
    return startTime + chrono::milliseconds(1700) < chrono::high_resolution_clock::now();
#endif
}

int mod(int a, int m) {
    return (a+m)%m;
}

template <typename T>
void printVector(vector<T> v) {
    for (auto i: v) {
        cerr << i << " ";
    }
    cerr << endl;
}

class TravelingSalesmanProblem{
public:
    int points;
    double *x;
    double *y;
    double **distance;
    int **neighbors;
    
    template <class T>
    static TravelingSalesmanProblem *createFromStream(T &cin){
        TravelingSalesmanProblem *instance = new TravelingSalesmanProblem();
        cin >> instance->points;
        instance->y = new double[instance->points];
        instance->x = new double[instance->points];
        for (auto i = 0; i < instance->points; i++) {
            cin >> instance->x[i] >> instance->y[i];
        }
        initInstance(instance);
        return instance;
    }
    
    static TravelingSalesmanProblem *createFromStdin(){
        return createFromStream(cin);
    }
    
    static TravelingSalesmanProblem *createFromFile() {
        fstream file("tsp.txt");
        return createFromStream(file);
    }
    
    static TravelingSalesmanProblem *createRandom(int size) {
        TravelingSalesmanProblem *instance = new TravelingSalesmanProblem();
        instance->points = size;
        instance->y = new double[instance->points];
        instance->x = new double[instance->points];
        for (auto i = 0; i < instance->points; i++) {
            instance->x[i] = rand() % 1000;
            instance->y[i] = rand() % 1000;
        }
        initInstance(instance);
        return instance;
    }
    
    static void initInstance(TravelingSalesmanProblem *instance) {
        instance->distance = new double*[instance->points];
        instance->neighbors = new int*[instance->points];
        
        for (int i = 0; i < instance->points; i++) {
            instance->distance[i] = new double[instance->points];
            for (int j = 0; j < instance->points; j++) {
                instance->distance[i][j] = sqrt((instance->x[i]-instance->x[j])*(instance->x[i]-instance->x[j]) + (instance->y[i]-instance->y[j])*(instance->y[i]-instance->y[j]));
            }
        }
        
        for (int i = 0; i < instance->points; i++) {
            instance->neighbors[i] = new int[instance->points];
            for (int j = 0; j < instance->points; j++) {
                instance->neighbors[i][j] = j;
            }
            sort(instance->neighbors[i], instance->neighbors[i]+instance->points, [instance, i](int a, int b){
                return instance->dist(i, a) < instance->dist(i, b);
            });
            assert(instance->dist(0, instance->neighbors[0][0]) <= instance->dist(0, instance->neighbors[0][instance->points-1]));
            
            //            for (int j = 1; j < instance->points; j++) {
            //                int k = j-1;
            //                int jValue = instance->neighbors[i][j];
            //                int jDistance = instance->dist(i, jValue);
            //                while (k >= 0 && instance->dist(i, instance->neighbors[i][k]) > jDistance) {
            //                    instance->neighbors[i][k+1] = instance->neighbors[i][k];
            //                    k--;
            //                }
            //                instance->neighbors[i][k+1] = j;
            //            }
        }
        
        //                for (int i = 0; i < instance->points; i++) {
        //                    cerr << "Neighbors of " << i << endl;
        //                    for (int j = 0; j < instance->points; j++) {
        //                        int neighbor = instance->neighbors[i][j];
        //                        cerr << "(" << neighbor << ", " << instance->dist(i, neighbor) << ") " ;
        //                    }
        //                    cerr << endl;
        //                }
        
    }
    
    static TravelingSalesmanProblem *testInstance(){
        TravelingSalesmanProblem *instance = new TravelingSalesmanProblem();
        instance->points = 10;
        instance->x = new double[instance->points];
        instance->y = new double[instance->points];
        instance->x[0] = 95.0129;        //    95.0129 61.5432
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
        initInstance(instance);
        return instance;
    }
    
    inline double dist(int a, int b){
        return distance[a][b];
    }
    
    void printBits(unsigned long bits) {
        for (int i = 31; i >= 0; i--) {
            cout << ((bits >> i) & 1) << " ";
        }
        cout << endl;
    }
    
    vector<int> nearestNeighbor(int start) {
        vector<int> vector;
        int *tour = new int[points];
        bool *used = new bool[points];
        for (int i = 0; i < points; i++) {
            used[i] = false;
        }
        
        tour[0] = start;
        used[start] = true;
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
            //            cout << tour[i] << endl;
        }
        return vector;
    }
    
    vector<int> greedy() {
        vector<int> tour(points, 0);
        vector<int> next(points, -1);
        vector<int> prev(points, -1);
        for (int i = 0; i < points; i++) {
            int best = -1;
            double bestValue = 99999999;
            int bestN = -1;
            for (int j = 0; j < points; j++) {
                if (next[j] != -1) {
                    continue;
                }
                for (int n2 = 1; n2 < points; n2++) {
                    auto n = neighbors[j][n2];
                    if (n == j) {
                        continue;
                    }
                    if (prev[n] != -1) {
                        continue;
                    }
                    
                    auto p = j;
                    
                    bool isCyclic = false;
                    while (prev[p] != -1) {
                        p = prev[p];
                        if (p == n) {
                            isCyclic = true;
                            break;
                        }
                    }
                    
                    
                    if (isCyclic && i != tour.size()-1) {
                        continue;
                    }
                    auto distance = dist(j, n);
                    if (distance < bestValue) {
                        best = j;
                        bestValue = distance;
                        bestN = n;
                    }
                    break;
                }
            }
            auto neighbor = bestN;
            next[best] = neighbor;
            prev[neighbor] = best;
        }
        for (int i = 0; i < points; i++) {
            assert(find(next.begin(), next.end(), i) != next.end());
        }
        
        int path = 0;
        for (int i = 0; i < points; i++) {
            tour[i] = path;
            path = next[path];
        }
        
        return tour;
    }
    
    vector<int> clarkeWright() {
        vector<int> tour(points, 0);
        vector<int> next(points, -1);
        vector<int> prev(points, -1);
        
        int hub = rand() % tour.size();
        
        for (int i = 0; i < points; i++) {
            int best = -1;
            double bestValue = -1;
            int bestN = -1;
            for (int j = 0; j < points; j++) {
                if (next[j] != -1) {
                    continue;
                }
                for (int n2 = 1; n2 < points; n2++) {
                    auto n = neighbors[j][n2];
                    if (n == j) {
                        continue;
                    }
                    if (prev[n] != -1) {
                        continue;
                    }
                    auto p = j;
                    bool isCyclic = false;
                    while (prev[p] != -1) {
                        p = prev[p];
                        if (p == n) {
                            isCyclic = true;
                            break;
                        }
                    }
                    
                    if (isCyclic && i != tour.size()-1) {
                        continue;
                    }
                    auto distance = (dist(hub, n) + dist(hub, j)) - (dist(j, n));
                    if (distance > bestValue) {
                        best = j;
                        bestValue = distance;
                        bestN = n;
                    }
                }
            }
            auto neighbor = bestN;
            next[best] = neighbor;
            prev[neighbor] = best;
        }
        for (int i = 0; i < points; i++) {
            assert(find(next.begin(), next.end(), i) != next.end());
        }
        
        int path = 0;
        for (int i = 0; i < points; i++) {
            tour[i] = path;
            path = next[path];
        }
        
        return tour;
    }
    
    
    double tourDistance(const vector<int> &tour) {
        double distance = 0;
        for (int i = 1; i < tour.size(); i++) {
            distance += dist(tour[i-1], tour[i]);
        }
        return distance + dist(tour[0], tour[tour.size()-1]);
    }
    
    void kopt2(vector<int> &tour) {
        bool didImprove = true;
        while (didImprove) {
            didImprove = false;
            for (int i = 0; i < tour.size(); i++) {
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return;
                }
                for (int j = i+2; j < tour.size(); j++) {
                    int nextJ = j+1 < tour.size() ? j+1 : 0;
                    double prevDistance = dist(tour[i], tour[i+1]) + dist(tour[j], tour[nextJ]);
                    double afterDistance = dist(tour[i], tour[j]) + dist(tour[i+1], tour[nextJ]);
                    if (afterDistance < prevDistance) {
                        reverse(&tour[i+1], &tour[j]+1);
                        didImprove = true;
                    }
                }
            }
        }
        return;
    }
    
    bool kopt2(vector<int> &tour, int i) {
        bool didImprove = false;
        for (int j = i+2; j < tour.size(); j++) {
            int nextJ = j+1 < tour.size() ? j+1 : 0;
            double prevDistance = dist(tour[i], tour[i+1]) + dist(tour[j], tour[nextJ]);
            double afterDistance = dist(tour[i], tour[j]) + dist(tour[i+1], tour[nextJ]);
            if (afterDistance < prevDistance) {
                didImprove = true;
                reverse(&tour[i+1], &tour[j]+1);
            }
        }
        return didImprove;
    }
    
    
    void kopt2neighbors(vector<int> &tour) {
        vector<int> reverseTour(tour.size(), 0);
        for (auto i = 0; i < tour.size(); i++) {
            reverseTour[tour[i]] = i;
        }
        while (true) {
            bool didSwap = false;
            for (int i = 0; i < tour.size(); i++) {
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return;
                }
                for (int n = 1; n < min(1000, points); n++) {
                    int j = reverseTour[neighbors[tour[i]][n]];
                    //                    if (dist(tour[i], tour[j]) > dist(tour[i], tour[i+1]) + 10) {
                    //                        break;
                    //                    }
                    if (j >= i+2) {
                        int nextJ = j+1 < tour.size() ? j+1 : 0;
                        double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[nextJ]);
                        double afterDistance2 = dist(tour[i], tour[j]) + dist(tour[i+1], tour[nextJ]);
                        if (afterDistance2 + 0.1 < prevDistance2) {
                            reverse(&tour[i]+1, &tour[j]+1);
                            for (auto k = i+1; k < j+1; k++) {
                                reverseTour[tour[k]] = k;
                            }
                            didSwap = true;
                        }
                    }
                }
                
            }
            if (!didSwap) {
                break;
            }
        }
        return;
    }
    
    void kopt2neighborsopt(vector<int> &tour) {
        vector<int> reverseTour(tour.size(), 0);
        for (auto i = 0; i < tour.size(); i++) {
            reverseTour[tour[i]] = i;
        }
        while (true) {
            bool didSwap = false;
            for (int i = 0; i < tour.size()-1; i++) {
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return;
                }
                for (int n = 1; n < min(100, points); n++) {
                    int j = reverseTour[neighbors[tour[i]][n]];
                    if (dist(tour[i], tour[i + 1]) < dist(tour[i], tour[j]))  break;
                    if (j >= i + 2 || j <= i - 2) {
                        int nextJ = j+1 < tour.size() ? j+1 : 0;
                        double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[nextJ]);
                        double afterDistance2 = dist(tour[i], tour[j]) + dist(tour[i+1], tour[nextJ]);
                        if (afterDistance2 + 0.1 < prevDistance2) {
                            int i2 = i;
                            int j2 = j;
                            if (i2 > j2)
                                swap(i2, j2);
                            reverse(&tour[i2] + 1, &tour[j2] + 1);
                            for (auto k = i2 + 1; k < j2 + 1; k++) {
                                reverseTour[tour[k]] = k;
                            }
                            didSwap = true;
                        }
                    }
                }
                
            }
            if (!didSwap) {
                break;
            }
        }
        return;
    }
    
    
    bool kopt2neighbors(vector<int> &tour, vector<int> &reverseTour, const int i) {
        bool didImprove = false;
        for (int n = 1; n < min(10, points); n++) {
            int jValue = neighbors[tour[i]][n];
            int j = reverseTour[jValue];
            if (j >= i+2) {
                int nextJ = j+1 < tour.size() ? j+1 : 0;
                //                        if (tour[i]+tour[i+1] < tour[i] + tour[j]) {
                //                            break;
                //                        }
                double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[nextJ]);
                double afterDistance2 = dist(tour[i], tour[j]) + dist(tour[i+1], tour[nextJ]);
                if (afterDistance2 < prevDistance2) {
                    reverse(&tour[i+1], &tour[j]+1);
                    for (auto k = i+1; k < j+1; k++) {
                        reverseTour[tour[k]] = k;
                    }
                    didImprove = true;
                }
            }
        }
        return didImprove;
    }
    
    
    void kopt3neighbors3(vector<int> &tour, int neighborLimit) {
        vector<int> reverseTour(tour.size(), 0);
        for (auto x = 0; x < tour.size(); x++) {
            reverseTour[tour[x]] = x;
        }
        
        bool didImprove = true;
        double distanceThreshold = 0.1;
        while (didImprove) {
            didImprove = false;
            for (int i = 0; i < tour.size()-1; i++) {
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return;
                }
                // First move
                for (int iNeighbor = 1; iNeighbor < min(neighborLimit, (int)tour.size()); iNeighbor++) {
                    int j = reverseTour[neighbors[tour[i]][iNeighbor]];
                    for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                        if (j <= i) {
                            break;
                        }
                        int k = reverseTour[neighbors[tour[i+1]][jNeighbor]];
                        if (j < k) {
                            assert(i < j && j < k);
                            assert(i >= 0);
                            assert(k < tour.size());
                            if (timeIsRunningOut()) {
                                cerr << "Time out" << endl;
                                return;
                            }
                            int nextK = k+1 < tour.size() ? k+1 : 0;
                            double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[j+1]) + dist(tour[k], tour[nextK]);
                            double afterDistance2 = dist(tour[i], tour[j]) + dist(tour[i+1], tour[k]) + dist(tour[j+1], tour[nextK]);
                            if (afterDistance2 + distanceThreshold < prevDistance2) {
                                auto bef = tourDistance(tour);
                                reverse(&tour[i]+1, &tour[j]+1);
                                reverse(&tour[j]+1, &tour[k]+1);
                                auto af = tourDistance(tour);
                                assert(bef > af);
                                didImprove = true;
                                for (auto x = 0; x < tour.size(); x++) {
                                    reverseTour[tour[x]] = x;
                                }
                            }
                        }
                    }
                    if (timeIsRunningOut()) {
                        cerr << "Time out" << endl;
                        return;
                    }
                    
                    // Second move
                    j = reverseTour[neighbors[tour[i]][iNeighbor]] - 1;
                    for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                        if (j <= i) {
                            break;
                        }
                        int k = reverseTour[neighbors[tour[j]][jNeighbor]] - 1;
                        if (k == -1) {
                            k = (int)tour.size()-1;
                        }
                        if (j < k) {
                            assert(i < j && j < k);
                            int nextK = k+1 < tour.size() ? k+1 : 0;
                            double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[j+1]) + dist(tour[k], tour[nextK]);
                            double afterDistance2 = dist(tour[i], tour[j+1]) + dist(tour[j], tour[nextK]) + dist(tour[k], tour[i+1]);
                            if (afterDistance2 + distanceThreshold < prevDistance2) {
                                auto bef = tourDistance(tour);
                                vector<int> tmp(tour.size(), -1);
                                int index = 0;
                                for (int x = 0; x <= i; x++) {
                                    tmp[x] = tour[x];
                                    index++;
                                }
                                for (int x = j+1; x <= k; x++) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                for (int x = i+1; x <= j; x++) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                for (int x = k+1; x < tour.size(); x++) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                
                                tour = tmp;
                                didImprove = true;
                                auto af = tourDistance(tour);
                                //                                assert(bef > af);
                                for (auto x = 0; x < tour.size(); x++) {
                                    reverseTour[tour[x]] = x;
                                }
                            }
                        }
                    }
                    if (timeIsRunningOut()) {
                        cerr << "Time out" << endl;
                        return;
                    }
                    
                    
                    // Third move
                    j = reverseTour[neighbors[tour[i+1]][iNeighbor]] - 1;
                    for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                        if (j <= i) {
                            break;
                        }
                        int k = reverseTour[neighbors[tour[i]][jNeighbor]];
                        if (k == -1) {
                            k = (int)tour.size()-1;
                        }
                        if (j < k) {
                            assert(i < j && j < k);
                            int nextK = k+1 < tour.size() ? k+1 : 0;
                            double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[j+1]) + dist(tour[k], tour[nextK]);
                            double afterDistance2 = dist(tour[i], tour[k]) + dist(tour[j+1], tour[i+1]) + dist(tour[j], tour[nextK]);
                            if (afterDistance2 + distanceThreshold< prevDistance2) {
                                auto bef = tourDistance(tour);
                                vector<int> tmp(tour.size());
                                int index = 0;
                                for (int x = 0; x <= i; x++) {
                                    tmp[x] = tour[x];
                                    index++;
                                }
                                for (int x = k; x >= j+1; x--) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                for (int x = i+1; x <= j; x++) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                for (int x = k+1; x < tour.size(); x++) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                tour = tmp;
                                assert(bef > tourDistance(tour));
                                didImprove = true;
                                for (auto x = 0; x < tour.size(); x++) {
                                    reverseTour[tour[x]] = x;
                                }
                            }
                        }
                        if (timeIsRunningOut()) {
                            cerr << "Time out" << endl;
                            return;
                        }
                        
                        
                        // Fourth move
                        j = reverseTour[neighbors[tour[i]][iNeighbor]] - 1;
                        for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                            if (j <= i) {
                                break;
                            }
                            int k = reverseTour[neighbors[tour[j]][jNeighbor]];
                            if (k == -1) {
                                k = (int)tour.size()-1;
                            }
                            if (j < k) {
                                assert(i < j && j < k);
                                int nextK = k+1 < tour.size() ? k+1 : 0;
                                double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[j+1]) + dist(tour[k], tour[nextK]);
                                double afterDistance2 = dist(tour[i], tour[j+1]) + dist(tour[k], tour[j]) + dist(tour[i+1], tour[nextK]);
                                if (afterDistance2 + distanceThreshold < prevDistance2) {
                                    auto bef = tourDistance(tour);
                                    vector<int> tmp(tour.size());
                                    int index = 0;
                                    for (int x = 0; x <= i; x++) {
                                        tmp[x] = tour[x];
                                        index++;
                                    }
                                    for (int x = j+1; x <= k; x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    for (int x = j; x >= i+1; x--) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    for (int x = k+1; x < tour.size(); x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    tour = tmp;
                                    assert(bef > tourDistance(tour));
                                    didImprove = true;
                                    for (auto x = 0; x < tour.size(); x++) {
                                        reverseTour[tour[x]] = x;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    void kopt3neighbors2(vector<int> &tour, int neighborLimit) {
        vector<int> reverseTour(tour.size(), 0);
        for (auto x = 0; x < tour.size(); x++) {
            reverseTour[tour[x]] = x;
        }
        
        bool didImprove = true;
        double distanceThreshold = 0.1;
        while (didImprove) {
            didImprove = false;
            for (int i = 0; i < tour.size()-1; i++) {
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return;
                }
                // First move
                for (int iNeighbor = 1; iNeighbor < min(neighborLimit, (int)tour.size()); iNeighbor++) {
                    int j = reverseTour[neighbors[tour[i]][iNeighbor]];
                    if (j <= i) {
                        continue;
                    }
                    for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                        int k = reverseTour[neighbors[tour[i+1]][jNeighbor]];
                        if (j < k) {
                            assert(i < j && j < k);
                            assert(i >= 0);
                            assert(k < tour.size());
                            if (timeIsRunningOut()) {
                                cerr << "Time out" << endl;
                                return;
                            }
                            int nextK = k+1 < tour.size() ? k+1 : 0;
                            double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[j+1]) + dist(tour[k], tour[nextK]);
                            double afterDistance2 = dist(tour[i], tour[j]) + dist(tour[i+1], tour[k]) + dist(tour[j+1], tour[nextK]);
                            if (afterDistance2 + distanceThreshold < prevDistance2) {
                                auto bef = tourDistance(tour);
                                reverse(&tour[i]+1, &tour[j]+1);
                                reverse(&tour[j]+1, &tour[k]+1);
                                auto af = tourDistance(tour);
                                assert(bef > af);
                                didImprove = true;
                                for (auto x = 0; x < tour.size(); x++) {
                                    reverseTour[tour[x]] = x;
                                }
                            }
                        }
                    }
                }
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return;
                }
                
                // Second move
                for (int iNeighbor = 1; iNeighbor < min(neighborLimit, (int)tour.size()); iNeighbor++) {
                    int j = reverseTour[neighbors[tour[i]][iNeighbor]] - 1;
                    if (j <= i) {
                        continue;
                    }
                    for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                        int k = reverseTour[neighbors[tour[j]][jNeighbor]] - 1;
                        if (k == -1) {
                            k = (int)tour.size()-1;
                        }
                        if (j < k) {
                            assert(i < j && j < k);
                            int nextK = k+1 < tour.size() ? k+1 : 0;
                            double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[j+1]) + dist(tour[k], tour[nextK]);
                            double afterDistance2 = dist(tour[i], tour[j+1]) + dist(tour[j], tour[nextK]) + dist(tour[k], tour[i+1]);
                            if (afterDistance2 + distanceThreshold < prevDistance2) {
                                auto bef = tourDistance(tour);
                                vector<int> tmp(tour.size(), -1);
                                int index = 0;
                                for (int x = 0; x <= i; x++) {
                                    tmp[x] = tour[x];
                                    index++;
                                }
                                for (int x = j+1; x <= k; x++) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                for (int x = i+1; x <= j; x++) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                for (int x = k+1; x < tour.size(); x++) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                
                                tour = tmp;
                                didImprove = true;
                                auto af = tourDistance(tour);
                                //                                assert(bef > af);
                                for (auto x = 0; x < tour.size(); x++) {
                                    reverseTour[tour[x]] = x;
                                }
                            }
                        }
                    }
                }
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return;
                }
                
                
                // Third move
                for (int iNeighbor = 1; iNeighbor < min(neighborLimit, (int)tour.size()); iNeighbor++) {
                    int j = reverseTour[neighbors[tour[i+1]][iNeighbor]] - 1;
                    if (j <= i) {
                        continue;
                    }
                    for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                        int k = reverseTour[neighbors[tour[i]][jNeighbor]];
                        if (k == -1) {
                            k = (int)tour.size()-1;
                        }
                        if (j < k) {
                            assert(i < j && j < k);
                            int nextK = k+1 < tour.size() ? k+1 : 0;
                            double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[j+1]) + dist(tour[k], tour[nextK]);
                            double afterDistance2 = dist(tour[i], tour[k]) + dist(tour[j+1], tour[i+1]) + dist(tour[j], tour[nextK]);
                            if (afterDistance2 + distanceThreshold< prevDistance2) {
                                auto bef = tourDistance(tour);
                                vector<int> tmp(tour.size());
                                int index = 0;
                                for (int x = 0; x <= i; x++) {
                                    tmp[x] = tour[x];
                                    index++;
                                }
                                for (int x = k; x >= j+1; x--) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                for (int x = i+1; x <= j; x++) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                for (int x = k+1; x < tour.size(); x++) {
                                    tmp[index] = tour[x];
                                    index++;
                                }
                                tour = tmp;
                                assert(bef > tourDistance(tour));
                                didImprove = true;
                                for (auto x = 0; x < tour.size(); x++) {
                                    reverseTour[tour[x]] = x;
                                }
                            }
                        }
                    }
                    if (timeIsRunningOut()) {
                        cerr << "Time out" << endl;
                        return;
                    }
                    
                    // Fourth move
                    for (int iNeighbor = 1; iNeighbor < min(neighborLimit, (int)tour.size()); iNeighbor++) {
                        int j = reverseTour[neighbors[tour[i]][iNeighbor]] - 1;
                        if (j <= i) {
                            continue;
                        }
                        for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                            int k = reverseTour[neighbors[tour[j]][jNeighbor]];
                            if (k == -1) {
                                k = (int)tour.size()-1;
                            }
                            if (j < k) {
                                assert(i < j && j < k);
                                int nextK = k+1 < tour.size() ? k+1 : 0;
                                double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[j+1]) + dist(tour[k], tour[nextK]);
                                double afterDistance2 = dist(tour[i], tour[j+1]) + dist(tour[k], tour[j]) + dist(tour[i+1], tour[nextK]);
                                if (afterDistance2 + distanceThreshold < prevDistance2) {
                                    auto bef = tourDistance(tour);
                                    vector<int> tmp(tour.size());
                                    int index = 0;
                                    for (int x = 0; x <= i; x++) {
                                        tmp[x] = tour[x];
                                        index++;
                                    }
                                    for (int x = j+1; x <= k; x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    for (int x = j; x >= i+1; x--) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    for (int x = k+1; x < tour.size(); x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    tour = tmp;
                                    assert(bef > tourDistance(tour));
                                    didImprove = true;
                                    for (auto x = 0; x < tour.size(); x++) {
                                        reverseTour[tour[x]] = x;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    void kopt3neighbors2opt(vector<int> &tour, int neighborLimit) {
        vector<int> reverseTour(tour.size(), 0);
        for (auto x = 0; x < tour.size(); x++) {
            reverseTour[tour[x]] = x;
        }
        
        bool didImprove = true;
        double distanceThreshold = 0.1;
        while (didImprove) {
            didImprove = false;
            for (int i = 0; i < tour.size()-1; i++) {
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return;
                }
                // First move
                for (int iNeighbor = 1; iNeighbor < min(neighborLimit, (int)tour.size()); iNeighbor++) {
                    int j = reverseTour[neighbors[tour[i]][iNeighbor]];
                    if (dist(tour[i], tour[i + 1]) < dist(tour[i], tour[j])) break;
                    if (j >= i + 2 || j <= i - 2) {
                        for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                            int k = reverseTour[neighbors[tour[i+1]][jNeighbor]];
                            if (j < k && i < k) {
                                int nextK = k+1 < tour.size() ? k+1 : 0;
                                int i2 = i;
                                int j2 = j;
                                if (i2 > j2)
                                    swap(i2, j2);
                                double prevDistance2 = dist(tour[i2], tour[i2+1]) + dist(tour[j2], tour[j2+1]) + dist(tour[k], tour[nextK]);
                                double afterDistance2 = dist(tour[i2], tour[j2]) + dist(tour[i2+1], tour[k]) + dist(tour[j2+1], tour[nextK]);
                                if (afterDistance2 + distanceThreshold < prevDistance2) {
                                    auto bef = tourDistance(tour);
                                    reverse(&tour[i2]+1, &tour[j2]+1);
                                    reverse(&tour[j2]+1, &tour[k]+1);
                                    auto af = tourDistance(tour);
                                    assert(bef > af);
                                    didImprove = true;
                                    for (auto x = 0; x < tour.size(); x++) {
                                        reverseTour[tour[x]] = x;
                                    }
                                }
                            }
                        }
                    }
                }
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return;
                }
                
                
                // Second move
                for (int iNeighbor = 1; iNeighbor < min(neighborLimit, (int)tour.size()); iNeighbor++) {
                    int j = reverseTour[neighbors[tour[i]][iNeighbor]] - 1;
                    if (j == -1) continue;
                    if (dist(tour[i], tour[i + 1]) < dist(tour[i], tour[j+1])) break;
                    if (j >= i + 2 || j <= i - 2) {
                        for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                            int k = reverseTour[neighbors[tour[j]][jNeighbor]] - 1;
                            if (k == -1) {
                                k = (int)tour.size()-1;
                            }
                            if (j < k && i < k) {
                                int nextK = k+1 < tour.size() ? k+1 : 0;
                                int i2 = i;
                                int j2 = j;
                                if (i2 > j2)
                                    swap(i2, j2);
                                double prevDistance2 = dist(tour[i2], tour[i2+1]) + dist(tour[j2], tour[j2+1]) + dist(tour[k], tour[nextK]);
                                double afterDistance2 = dist(tour[i2], tour[j2+1]) + dist(tour[j2], tour[nextK]) + dist(tour[k], tour[i2+1]);
                                if (afterDistance2 + distanceThreshold < prevDistance2) {
                                    auto bef = tourDistance(tour);
                                    vector<int> tmp(tour.size(), -1);
                                    int index = 0;
                                    for (int x = 0; x <= i2; x++) {
                                        tmp[x] = tour[x];
                                        index++;
                                    }
                                    for (int x = j2+1; x <= k; x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    for (int x = i2+1; x <= j2; x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    for (int x = k+1; x < tour.size(); x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    
                                    tour = tmp;
                                    didImprove = true;
                                    auto af = tourDistance(tour);
                                    //                                assert(bef > af);
                                    for (auto x = 0; x < tour.size(); x++) {
                                        reverseTour[tour[x]] = x;
                                    }
                                }
                            }
                        }
                    }
                }
                
                // Third move
                for (int iNeighbor = 1; iNeighbor < min(neighborLimit, (int)tour.size()); iNeighbor++) {
                    int j = reverseTour[neighbors[tour[i+1]][iNeighbor]] - 1;
                    if (j == -1) continue;
                    if (dist(tour[i], tour[i + 1]) < dist(tour[i+1], tour[j+1])) break;
                    if (j >= i + 2 || j <= i - 2) {
                        for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                            int k = reverseTour[neighbors[tour[i]][jNeighbor]];
                            if (k == -1) {
                                k = (int)tour.size()-1;
                            }
                            if (j < k && i < k) {
                                
                                int nextK = k+1 < tour.size() ? k+1 : 0;
                                int i2 = i;
                                int j2 = j;
                                if (i2 > j2)
                                    swap(i2, j2);
                                double prevDistance2 = dist(tour[i2], tour[i2+1]) + dist(tour[j2], tour[j2+1]) + dist(tour[k], tour[nextK]);
                                double afterDistance2 = dist(tour[i2], tour[k]) + dist(tour[j2+1], tour[i2+1]) + dist(tour[j2], tour[nextK]);
                                if (afterDistance2 + distanceThreshold< prevDistance2) {
                                    auto bef = tourDistance(tour);
                                    vector<int> tmp(tour.size());
                                    int index = 0;
                                    for (int x = 0; x <= i2; x++) {
                                        tmp[x] = tour[x];
                                        index++;
                                    }
                                    for (int x = k; x >= j2+1; x--) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    for (int x = i2+1; x <= j2; x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    for (int x = k+1; x < tour.size(); x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    tour = tmp;
                                    assert(bef > tourDistance(tour));
                                    didImprove = true;
                                    for (auto x = 0; x < tour.size(); x++) {
                                        reverseTour[tour[x]] = x;
                                    }
                                }
                            }
                        }
                    }
                }
                
                // Fourth move
                for (int iNeighbor = 1; iNeighbor < min(neighborLimit, (int)tour.size()); iNeighbor++) {
                    int j = reverseTour[neighbors[tour[i]][iNeighbor]] - 1;
                    if (j == -1) continue;
                    if (dist(tour[i], tour[i + 1]) < dist(tour[i], tour[j+1])) break;
                    if (j >= i + 2 || j <= i - 2) {
                        for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                            int k = reverseTour[neighbors[tour[j]][jNeighbor]];
                            if (k == -1) {
                                k = (int)tour.size()-1;
                            }
                            if (j < k && i < k) {
                                int nextK = k+1 < tour.size() ? k+1 : 0;
                                int i2 = i;
                                int j2 = j;
                                if (i2 > j2)
                                    swap(i2, j2);
                                double prevDistance2 = dist(tour[i2], tour[i2+1]) + dist(tour[j2], tour[j2+1]) + dist(tour[k], tour[nextK]);
                                double afterDistance2 = dist(tour[i2], tour[j2+1]) + dist(tour[k], tour[j2]) + dist(tour[i2+1], tour[nextK]);
                                if (afterDistance2 + distanceThreshold < prevDistance2) {
                                    auto bef = tourDistance(tour);
                                    vector<int> tmp(tour.size());
                                    int index = 0;
                                    for (int x = 0; x <= i2; x++) {
                                        tmp[x] = tour[x];
                                        index++;
                                    }
                                    for (int x = j2+1; x <= k; x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    for (int x = j2; x >= i2+1; x--) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    for (int x = k+1; x < tour.size(); x++) {
                                        tmp[index] = tour[x];
                                        index++;
                                    }
                                    tour = tmp;
                                    assert(bef > tourDistance(tour));
                                    didImprove = true;
                                    for (auto x = 0; x < tour.size(); x++) {
                                        reverseTour[tour[x]] = x;
                                    }
                                }
                            }
                        }
                    }
                }
                
            }
        }
    }
    
    
    vector<int> linKernighanMove3(const vector<int> &move, int moveSize, const vector<int> &tour, const vector<int> &reverseTour) {
        vector<int> inverseSubtours(tour.size(), -1);
        vector<pair<int, int>> subtourIndexes;
        vector<int> splitIndexes;
        splitIndexes.reserve(10);
        
        for (int i = 0; i < moveSize; i += 2) {
            splitIndexes.push_back((move[i] + 1) % tour.size() == move[i + 1] ? move[i] : move[i + 1]);
        }
        if (splitIndexes[(int)splitIndexes.size() - 1] != tour.size() - 1) {
            splitIndexes.push_back((int)tour.size() - 1);
        }
        sort(splitIndexes.begin(), splitIndexes.end());
        
        int lastSplitIndex = -1;
        for (auto splitIndex : splitIndexes) {
            if (lastSplitIndex+1 <= splitIndex) {
                subtourIndexes.push_back(pair<int, int>(lastSplitIndex+1, splitIndex));
                inverseSubtours[lastSplitIndex + 1] = (int)subtourIndexes.size() - 1;
                inverseSubtours[splitIndex] = (int)subtourIndexes.size() - 1;
            }
            lastSplitIndex = splitIndex;
        }
        vector<bool> isSubtourUsed(subtourIndexes.size(), false);
        int lastIndex = subtourIndexes[0].second;
        isSubtourUsed[0] = true;
        int size = subtourIndexes[0].second - subtourIndexes[0].first + 1;
        while (size < tour.size()) {
            auto index = lastIndex;
            auto moveIndex = find(move.begin(), move.begin() + moveSize, index) - move.begin();
            if (moveIndex >= moveSize) return vector<int>();
            auto connectIndex = move[(moveIndex + (moveIndex % 2 == 0 ? moveSize - 1 : 1)) % moveSize];
            auto subtourIndex = inverseSubtours[connectIndex];
            if (subtourIndex == -1) {
                return vector<int>();
            }
            if (isSubtourUsed[subtourIndex]) return vector<int>();
            isSubtourUsed[subtourIndex] = true;
            if (subtourIndexes[subtourIndex].first == connectIndex) {
                lastIndex = subtourIndexes[subtourIndex].second;
            } else {
                lastIndex = subtourIndexes[subtourIndex].first;
            }
            size += subtourIndexes[subtourIndex].second - subtourIndexes[subtourIndex].first + 1;
        }
        
        for (int i = 0; i < isSubtourUsed.size(); i++) {
            isSubtourUsed[i] = false;
        }
        
        
        // Construct tour
        vector<int> result;
        result.reserve(tour.size());
        for (auto sub: subtourIndexes) {
            //            cout << sub.first << ", " << sub.second << endl;
        }
        result.insert(result.end(), tour.begin() + subtourIndexes[0].first, tour.begin() + subtourIndexes[0].second+1);
        isSubtourUsed[0] = true;
        
        while (result.size() < tour.size()) {
            auto index = reverseTour[result[result.size()-1]];
            auto moveIndex = find(move.begin(), move.begin() + moveSize, index) - move.begin();
            if (moveIndex >= moveSize) return vector<int>();
            auto connectIndex = move[(moveIndex + (moveIndex % 2 == 0 ? moveSize - 1 : 1)) % moveSize];
            auto subtourIndex = inverseSubtours[connectIndex];
            if (subtourIndex == -1) {
                return vector<int>();
            }
            if (isSubtourUsed[subtourIndex]) return vector<int>();
            isSubtourUsed[subtourIndex] = true;
            if (subtourIndexes[subtourIndex].first == connectIndex) {
                result.insert(result.end(), tour.begin() + subtourIndexes[subtourIndex].first, tour.begin() + subtourIndexes[subtourIndex].second + 1);
            } else {
                vector<int> tmp(tour.begin() + subtourIndexes[subtourIndex].first, tour.begin() + subtourIndexes[subtourIndex].second+1);
                result.insert(result.end(), tmp.rbegin(), tmp.rend());
            }
        }
        
        //        cout << "Result move" << endl;
        //        printVector(result);
        //        cout << "em" << endl;
        //        assert(result.size() == tour.size());
        //        for (auto isUsed : isSubtourUsed)
        //         {
        //         assert(isUsed);
        //         }
        
        return result;
    }
    
    
    
    // Best algorithm for tsp known to man
    void linKernighan2(vector<int> &tour) {
        vector<int> reverseTour(tour.size());
        for (int i = 0; i < tour.size(); i++) {
            reverseTour[tour[i]] = i;
        }
        int size = (int)tour.size();
        const auto next = [&](int i) { return mod(i + 1, size); };
        const auto prev = [&](int i) { return mod(i - 1, size); };
        bool didImprove = true;
        vector<int> t(size, -1);
        vector<int> tNeighbor(size, 0);
        vector<double> gains(size, 0);
        int i = 0;
        int iterations = 0;
        int moves = 0;
        int successfulMoves = 0;
        while (true) {
            iterations++;
                        if (timeIsRunningOut()) {
                            return;
                        }
            if (i % 2 == 0) {
                if (i > 8) {
                    i--;
                    continue;
                }
                if (i >= 2) {
                    tNeighbor[i]++;
                    if (tNeighbor[i] == 8) {
                        tNeighbor[i] = 0;
                        i--;
                        continue;
                    }
                    
                    int neighbor = this->neighbors[tour[t[i-1]]][tNeighbor[i]];
//                    t[i] = find(tour.begin(), tour.end(), neighbor) - tour.begin();
                    t[i] = reverseTour[neighbor];
                    if (t[i] == next(t[i - 1]) || t[i] == prev(t[i - 1]) || t[i] == t[i - 1]) {
                        continue;
                    }
                    
                    auto xDistance = dist(tour[t[i - 1]], tour[t[i - 2]]);
                    auto yDistance = dist(tour[t[i]], tour[t[i - 1]]);
                    auto gain = xDistance - yDistance;
                    gains[i] = gains[i-2] + gain;
                    if (gains[i] <= 0) {
                        tNeighbor[i] = 0;
                        i--;
                        continue;
                    }
                    
                } else {
                    t[i]++;
                    if (t[i] == tour.size()) {
                        break;
                    }
                }
            } else {
                if (t[i] == -1) {
                    t[i] = next(t[i - 1]);
                }
                else if (t[i] == next(t[i - 1])) {
                    t[i] = prev(t[i - 1]);
                }
                else {
                    t[i] = -1;
                    i--;
                    continue;
                }
                if (find(t.begin(), t.begin() + i, t[i]) != t.begin() + i) {
                    continue;
                }
                if (i >= 2) {
                    
                    // Ended tour??
                    if (t[i] == next(t[0]) || t[i] == prev(t[0]) || t[i] == t[0]) {
                        continue;
                    }
//                    double gain = 0;
//                    for (int j = 0; j < i; j++) {
//                        auto distance = dist(tour[t[j]], tour[t[j + 1]]);
//                        gain += distance * (j % 2 == 0 ? 1 : -1);
//                    }
                    
//                    double lastGain = -dist(tour[t[i]], tour[t[i-1]]);
//                    double gain2 = gain + lastGain;
                    
//                    double doubleThreshold = 0.00000001;
//                    assert(gain2 - gains[i-1] < doubleThreshold  && gain2 - gains[i-1] > -doubleThreshold);
//                    gain -= dist(tour[t[i]], tour[t[0]]);
                    
                    double tourDist = gains[i-1] + dist(tour[t[i]], tour[t[i-1]]) - dist(tour[t[i]], tour[t[0]]);
//                    assert(tourDist - gain < doubleThreshold && tourDist - gain > -doubleThreshold);
                    
                    // Check if better tour was constructed
                    if (tourDist > 0.000000001) {
                        auto potentialTour = linKernighanMove3(t, i+1, tour, reverseTour);
                        //                        auto potentialTour = linKernighanMove2(t, i+1, tour);
                        moves++;
                        if (potentialTour.size() != 0) {
                            auto bef = tourDistance(tour);
                            tour = potentialTour;
                            auto af = tourDistance(tour);
                            assert(bef >= af);
                            for (int i = 0; i < tour.size(); i++) {
                                reverseTour[tour[i]] = i;
                            }
                            didImprove = true;
                            for (int j = 0; j <= i; j++) {
                                t[j] = -1;
                                tNeighbor[j] = 0;
                            }
                            i = -1;
                            successfulMoves++;
                        } else {
                            //                            cerr << "Did not work" << endl;
                            //                            printVector(t);
                        }
                        
                    }
                }
            }
            i++;
            if (i == tour.size()) {
                i--;
            }
        }
        
        cerr << "Number of moves " << moves << endl;
        cerr << "Number of successful moves " << successfulMoves << endl;
        
    }
    
};


//int main(int argc, char **argv) {
//#ifdef DEBUG
//    auto instance = TravelingSalesmanProblem::createRandom(500);
//#else
//    auto instance = TravelingSalesmanProblem::createFromStdin();
//#endif
//    if (instance->points == 1) {
//        cout << 0;
//        return 0;
//    }
//    vector<int> greedyTour = instance->greedy();
//    random_shuffle(greedyTour.begin(), greedyTour.end());
//    vector<int> linKernighanTour = greedyTour;
//    auto kopt2Tour = linKernighanTour;
//    auto kopt3Tour = linKernighanTour;
//
//
//
//    instance->linKernighan2(linKernighanTour);
//    instance->kopt3neighbors2opt(linKernighanTour, 50);
//    instance->kopt2neighbors(linKernighanTour);
//
//
//
//    instance->kopt3neighbors2opt(kopt3Tour, 50);
//    instance->kopt2neighbors(kopt2Tour);
//
//
//    vector<int> minimumTour = linKernighanTour;
//
//    int iterations = 0;
//    double minimumDistance = 99999999999;
//
////        while (!timeIsRunningOut()) {
////            if (iterations > 10) {
////                break;
////            }
////
////            instance->kopt3neighbors2opt(greedyTour, 50);
////            instance->kopt2neighborsopt(greedyTour);
////            double distance = instance->tourDistance(greedyTour);
////            if (distance < minimumDistance) {
////                minimumDistance = distance;
////                minimumTour = greedyTour;
////            }
////
////            if (iterations < instance->points) {
////                greedyTour = instance->nearestNeighbor(iterations);
////            } else {
////                random_shuffle(greedyTour.begin(), greedyTour.end());
////            }
////            iterations++;
////        }
//
//    for (auto i: minimumTour) {
//        cout << i << endl;
//    }
//
//    for (int i = 0; i < instance->points; i++) {
//        assert(find(minimumTour.begin(), minimumTour.end(), i) != minimumTour.end());
//    }
//
//    cerr << "Greedy \t\t\t" << instance->tourDistance(greedyTour) << endl;
//    cerr << "2-opt\t\t\t" << instance->tourDistance(kopt2Tour) << endl;
//    cerr << "3-opt\t\t\t" << instance->tourDistance(kopt3Tour) << endl;
//    cerr << "Lin kernighan\t" << instance->tourDistance(linKernighanTour) << endl;
//    auto endTime = chrono::high_resolution_clock::now();
//    auto duration = (endTime - startTime);
//    cerr << "Time: " << duration.count() / 1e9 << endl;
////    cerr << "Iterations: " << iterations << endl;
//    delete instance;
//    return 0;
//}

int main(int argc, char **argv) {
#ifdef DEBUG
    auto instance = TravelingSalesmanProblem::createRandom(200);
#else
    auto instance = TravelingSalesmanProblem::createFromStdin();
#endif
    if (instance->points == 1) {
        cout << 0;
        return 0;
    }
    vector<int> greedyTour = instance->greedy();

    double minimumDistance = 99999999999;

//    if (instance->points < 400) {
        random_shuffle(greedyTour.begin(), greedyTour.end());
        instance->linKernighan2(greedyTour);
        minimumDistance = instance->tourDistance(greedyTour);
//    }
    
    vector<int> minimumTour = greedyTour;
    int iterations = 1;
    
//    while (!timeIsRunningOut()) {
//#ifdef DEBUG
//        if (iterations > 0) {
//            break;
//        }
//#endif
//        instance->kopt3neighbors2opt(greedyTour, 50);
//        instance->kopt2neighbors(greedyTour);
//        double distance = instance->tourDistance(greedyTour);
//        if (distance < minimumDistance) {
//            minimumDistance = distance;
//            minimumTour = greedyTour;
//        }
//        
//        if (iterations < instance->points) {
//            greedyTour = instance->nearestNeighbor(instance->points - iterations - 1);
//        } else {
//            random_shuffle(greedyTour.begin(), greedyTour.end());
//        }
//        iterations++;
//    }
//    
//    
    for (auto i: minimumTour) {
        cout << i << endl;
    }
    cerr << "Tour length: " << instance->tourDistance(minimumTour) << endl;
    auto endTime = chrono::high_resolution_clock::now();
    auto duration = (endTime - startTime);
    cerr << "Time: " << duration.count() / 1e9 << endl;
    delete instance;
    return 0;
}