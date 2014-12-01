#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <chrono>

using namespace std;
auto startTime = chrono::high_resolution_clock::now();

inline bool timeIsRunningOut() {
    return startTime + chrono::milliseconds(1500) < chrono::high_resolution_clock::now();
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
    
    static TravelingSalesmanProblem *createFromStdin(){
        TravelingSalesmanProblem *instance = new TravelingSalesmanProblem();;
        cin >> instance->points;
        instance->y = new double[instance->points];
        instance->x = new double[instance->points];
        for (auto i = 0; i < instance->points; i++)
            cin >> instance->x[i] >> instance->y[i];
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
            for (int j = 1; j < instance->points; j++) {
                int k = j-1;
                int jValue = instance->neighbors[i][j];
                int jDistance = instance->dist(i, jValue);
                while (instance->dist(i, instance->neighbors[i][k]) > jDistance && k >= 0) {
                    instance->neighbors[i][k+1] = instance->neighbors[i][k];
                    k--;
                }
                instance->neighbors[i][k+1] = j;
            }
        }
        
        //        for (int i = 0; i < instance->points; i++) {
        //            cerr << "Neighbors of " << i << endl;
        //            for (int j = 0; j < instance->points; j++) {
        //                int neighbor = instance->neighbors[i][j];
        //                cerr << "(" << neighbor << ", " << instance->dist(i, neighbor) << ") " ;
        //            }
        //            cerr << endl;
        //        }
        
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
            //            cout << tour[i] << endl;
        }
        return vector;
    }
    
    double tourDistance(vector<int> tour) {
        double distance = 0;
        for (int i = 1; i < tour.size(); i++) {
            distance += dist(tour[i-1], tour[i]);
        }
        return distance + dist(tour[0], tour[tour.size()-1]);
    }
    
    vector<int> kopt2(vector<int> &tour) {
        bool didImprove = true;
        while (didImprove) {
            didImprove = false;
            for (int i = 0; i < tour.size(); i++) {
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return tour;
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
        return tour;
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
    
    
    vector<int> kopt2neighbors(vector<int> tour) {
        vector<int> reverseTour(tour.size(), 0);
        for (auto i = 0; i < tour.size(); i++) {
            reverseTour[tour[i]] = i;
        }
        while (true) {
            bool didSwap = false;
            for (int i = 0; i < tour.size(); i++) {
                if (timeIsRunningOut()) {
                    cerr << "Time out" << endl;
                    return tour;
                }
                for (int n = 1; n < min(0, points); n++) {
                    int jValue = neighbors[i][n];
                    int j = reverseTour[jValue];
                    if (j >= i+2) {
                        int nextJ = j+1 < tour.size() ? j+1 : 0;
                        double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[nextJ]);
                        double afterDistance2 = dist(tour[i], tour[j]) + dist(tour[i+1], tour[nextJ]);
                        if (afterDistance2 < prevDistance2) {
                            reverse(&tour[i+1], &tour[j]+1);
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
        return tour;
    }
    
    
    vector<int> kopt3(vector<int> tour) {
        while (true) {
            bool didSwap = false;
            for (int i = 0; i < tour.size(); i++) {
                for (int j = i+2; j < tour.size(); j++) {
                    for (int k = j+2; k < tour.size(); k++) {
                        if (timeIsRunningOut()) {
                            cerr << "Time out" << endl;
                            return tour;
                        }
                        int nextK = k+1 < tour.size() ? k+1 : 0;
                        double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[j+1]) + dist(tour[k], tour[nextK]);
                        double afterDistance2 = dist(tour[i], tour[j]) + dist(tour[i+1], tour[k]) + dist(tour[j+1], tour[nextK]);
                        if (afterDistance2 < prevDistance2) {
                            reverse(&tour[i+1], &tour[j+1]);
                            reverse(&tour[j+1], &tour[nextK]);
                            didSwap = true;
                        }
                    }
                }
                
            }
            if (!didSwap) {
                break;
            }
        }
        return tour;
    }
};


int main(int argc, char **argv) {
#ifdef DEBUG
    auto instance = TravelingSalesmanProblem::testInstance();
#else
    auto instance = TravelingSalesmanProblem::createFromStdin();
#endif
    auto greedyTour = instance->greedy();
    auto tour = greedyTour;
    vector<int> minimumTour = tour;
    double minimumDistance = instance->tourDistance(tour);
    while (!timeIsRunningOut()) {
        tour = greedyTour;
        while(!timeIsRunningOut()) {
            if (!instance->kopt2(tour, rand() % instance->points)) {
                break;
            }
        }
        for (int i = 0; i < 6000; i++) {
            if (!timeIsRunningOut()) {
                instance->kopt2(tour, rand() % instance->points);
            }
        }
        double distance = instance->tourDistance(tour);
        if (distance < minimumDistance) {
            minimumDistance = distance;
            minimumTour = tour;
        }
//        random_shuffle(tour.begin(), tour.end());
    }
    
    cerr << instance->tourDistance(minimumTour) << endl;
    for (auto i: minimumTour) {
        cout << i << endl;
    }
    
    delete instance;
    return 0;
}