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
    return startTime + chrono::milliseconds(1800) < chrono::high_resolution_clock::now();
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
        instance->x[0] = 95.0129;        //	95.0129 61.5432
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
    
    
    vector<int> dynamicExact() {
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
        int node = 0;
        unsigned long S = (1 << n) - 1;
        for (int k = 1; k < n; k++) {
            double tmp = C[S][k] + dist(k, 0);
            opt = min(tmp, opt);
            node = k;
        }
        
        cout << "Opt: " << opt << endl;
        vector<int> tour;
        while (true) {
            tour.push_back(node);
            if (node == 0) {
                break;
            }
            int oldNode = node;
            node = parent[S][node];
            S = S^(1 << oldNode);
        }
        return tour;
    }
    
    // Best algorithm for tsp known to man
    void linKernighan(vector<int> &tour) {
        int size = (int) tour.size();
        const auto next = [&](int i) { return mod(i+1, size); };
        const auto prev = [&](int i) { return mod(i-1, size); };
        bool didImprove = true;
        
        vector<int> t(size, -1);
        for (t[0] = 0; t[0] < tour.size(); t[0]++) {
            didImprove = false;
            for (auto t1: {next(t[0]), prev(t[0])}) {
                t[1] = t1;
                int i = 2;
                for (t[i] = 0; t[i] < tour.size(); t[i]++) {
                    if (t[2] == next(t[1]) || t[2] == prev(t[1]) || t[2] == t[1]) { continue; }
                    int gain = dist(tour[t[0]], tour[t[1]]) - dist(tour[t[1]], tour[t[2]]);
                    if (gain < 0) { continue; };
                    for (auto t4: {next(t[i]), prev(t[i])}) {
                        t[i+1] = t4;
                        
                        // Ended tour??
                        if (t[i] == next(t[0]) || t[i] == prev(t[0])) { continue; }
                        int gain = 0;
                        for (int j = 0; j < i; j++) {
                            gain += dist(tour[t[j]], tour[t[j+1]]) * (j % 2 == 0 ? 1 : -1);
                        }
                        gain += dist(tour[t[i]], tour[t[0]]);
                        
                        // Check if better tour was constructed
                        if (gain > 0) {
                            auto bef = tourDistance(tour);
                            reverse(&tour[t[1]], &tour[t[3]]+1);
                            didImprove = true;
                            auto af = tourDistance(tour);
                            assert(bef >= af);
                        }
                    }
                    if (didImprove) {
                        t[0] = -1;
                        break;
                    }
                }
                if (didImprove) {
                    t[0] = -1;
                    break;
                }
            }
        }
    }
    
    //        while (true) {
    //            printVector(t);
    //            t[0]++;
    //            for (int i = 0; i < tour.size(); i++) {
    //                if (t[i] == tour.size()) {
    //                    t[i] = 0;
    //                    if (i+1 == tour.size()) {
    //                        goto hell;
    //                    }
    //                    t[i+1]++;
    //                } else {
    //                    break;
    //                }
    //            }
    //        }
    //        hell:
    //        cout << "Hell is done" << endl;
    
    
    vector<int> linKernighanMove(const vector<int> &move, const vector<int> &tour) {
        cout << "Moves" << endl;
        printVector(move);
        
        cout << "Tour" << endl;
        printVector(tour);
        vector<vector<int>> subtours;
        vector<int> inverseSubtours(tour.size(), -1);
        vector<int> splitIndexes;
        for (int i = 0; i < move.size(); i += 2) {
            splitIndexes.push_back((move[i]+1) % tour.size() == move[i+1] ? move[i] : move[i+1]);
        }
        if (splitIndexes[(int)splitIndexes.size()-1] != tour.size()-1) {
            splitIndexes.push_back((int)tour.size()-1);
        }
        sort(splitIndexes.begin(), splitIndexes.end());
        
        int lastSplitIndex = -1;
        for (auto splitIndex: splitIndexes) {
            vector<int> subtour(tour.begin() + (lastSplitIndex+1), tour.begin() + splitIndex + 1);
            lastSplitIndex = splitIndex;
            if (subtour.size() > 0) {
                subtours.push_back(subtour);
                inverseSubtours[lastSplitIndex+1] = (int)subtours.size() - 1;
                inverseSubtours[splitIndex] = (int)subtours.size() - 1;
            }
        }

        // Construct tour
        vector<int> result = subtours[0];
        while (result.size() < tour.size()) {
            auto index = find(tour.begin(), tour.end(), result[result.size()-1]) - tour.begin();
            auto moveIndex = find(move.begin(), move.end(), index) - move.begin();
            auto connectIndex = move[(moveIndex+(moveIndex%2==0 ? move.size()-1 : 1)) % move.size()];
            auto subtourIndex = inverseSubtours[connectIndex];
            if (subtourIndex == -1) {
                return result;
            }
            if (subtours[subtourIndex][0] == connectIndex) {
                result.insert(result.end(), subtours[subtourIndex].begin(), subtours[subtourIndex].end());
            } else {
                result.insert(result.end(), subtours[subtourIndex].rbegin(), subtours[subtourIndex].rend());
            }
        }
        cout << "Result move" << endl;
        printVector(result);
        
        return result;
    }
    
    
    // Best algorithm for tsp known to man
    void linKernighan3(vector<int> &tour) {
        int size = (int) tour.size();
        const auto next = [&](int i) { return mod(i+1, size); };
        const auto prev = [&](int i) { return mod(i-1, size); };
        bool didImprove = true;
        
        vector<int> t(size, -1);
        int i = 0;
        int iterations = 0;
        while (i < tour.size()) {
            
            iterations++;
            cout << "Iteration " << iterations << ":\t";
            printVector(t);
            if (i % 2 == 0) {
                if (i > 4) {
                    i--;
                    t[i]++;
                    continue;
                }
                
                t[i]++;
                if (t[i] == tour.size()) {
                    t[i] = -1;
                    if (i == 0) {
                        break;
                    }
                    i--;
                    continue;
                }
                if (i >= 2) {
                    int gain = dist(tour[t[i-1]], tour[t[i-2]]) - dist(tour[t[i]], tour[t[i-1]]);
                    if (t[i] == next(t[i-1]) || t[i] == prev(t[i-1]) || t[i] == t[i-1] || gain <= 0) {
                        continue;
                    } else {
                        cout << "Did find improvement" << endl;
                    }
                }
            } else {
                if (t[i] == -1) {
                    t[i] = next(t[i-1]);
                } else if (t[i] == next(t[i-1])) {
                    t[i] = prev(t[i-1]);
                } else {
                    t[i] = -1;
                    i--;
                    continue;
                }
                if (i >= 2) {
                    // Ended tour??
                    if (iterations == 239) {
                        
                    }
                    if (t[i] == next(t[0]) || t[i] == prev(t[0]) || t[i] == t[0]) {
                        cout << "Doesn't end tour" << endl;
                        continue;
                    }
                    int gain = 0;
                    for (int j = 0; j < i; j++) {
                        gain += dist(tour[t[j]], tour[t[j+1]]) * (j % 2 == 0 ? 1 : -1);
                    }
                    gain -= dist(tour[t[i]], tour[t[0]]);
                    
                    // Check if better tour was constructed
                    if (gain > 0) {
                        cout << "Did find improvement for entire tour " << i << endl;
                        vector<int> move;
                        for (int j = 0; j <= i; j++) {
                            move.push_back(t[j]);
                        }
                        
                        auto bef = tourDistance(tour);
                        
                        auto potentialTour = linKernighanMove(move, tour);
                        bool isTour = true;
                        for (int i = 0; i < tour.size(); i++) {
                            if (find(potentialTour.begin(), potentialTour.end(), i) == potentialTour.end()) {
                                isTour = false;
                            }
                        }
                        if (potentialTour.size() == tour.size() && isTour) {
                            cout << "DID MAKE A TOUR" << endl;
                            printVector(t);
                            printVector(tour);
                            tour = potentialTour;
                            didImprove = true;
                        }
                        
                        auto af = tourDistance(tour);
                        printVector(tour);
                        assert(bef >= af);
                    } else {
                        cout << "Did not gain from improvement" << endl;
                    }
                }
            }
            i++;
            if (i == tour.size()) {
                i--;
            }
        }
    }
    
};


int main(int argc, char **argv) {
#ifdef DEBUG
    auto instance = TravelingSalesmanProblem::createRandom(20);
#else
    auto instance = TravelingSalesmanProblem::createFromStdin();
#endif
    if (instance->points == 1) {
        cout << 0;
        return 0;
    }
    vector<int> greedyTour = instance->greedy();
    instance->linKernighan3(greedyTour);
    //            instance->kopt2neighborsopt(greedyTour);
    //
    vector<int> minimumTour = greedyTour;
    
    
    
    int iterations = 0;
    double minimumDistance = 99999999999;
    
    //    while (!timeIsRunningOut()) {
    //
    //        instance->kopt3neighbors2opt(greedyTour, 80);
    //        instance->kopt2neighborsopt(greedyTour);
    //        double distance = instance->tourDistance(greedyTour);
    //        if (distance < minimumDistance) {
    //            minimumDistance = distance;
    //            minimumTour = greedyTour;
    //        }
    //
    //        if (iterations < instance->points) {
    //            greedyTour = instance->nearestNeighbor(iterations);
    //        } else {
    //            random_shuffle(greedyTour.begin(), greedyTour.end());
    //        }
    //        iterations++;
    //    }
    //
    
    for (auto i: minimumTour) {
        cout << i << endl;
    }
    
    for (int i = 0; i < instance->points; i++) {
        assert(find(minimumTour.begin(), minimumTour.end(), i) != minimumTour.end());
    }
    
    
    
    
    cerr << "Tour length: " << instance->tourDistance(minimumTour) << endl;
    auto endTime = chrono::high_resolution_clock::now();
    auto duration = (endTime - startTime);
    cerr << "Time: " << duration.count() / 1e9 << endl;
    cerr << "Iterations: " << iterations << endl;
    delete instance;
    return 0;
}