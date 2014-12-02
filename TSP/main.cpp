#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <chrono>
#include <fstream>

using namespace std;
auto startTime = chrono::high_resolution_clock::now();

inline bool timeIsRunningOut() {
#ifdef DEBUG
    return false;
#else
    return startTime + chrono::milliseconds(1500) < chrono::high_resolution_clock::now();
#endif
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
        TravelingSalesmanProblem *instance = new TravelingSalesmanProblem();;
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
                while (k >= 0 && instance->dist(i, instance->neighbors[i][k]) > jDistance) {
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
                for (int n = 1; n < min(200, points); n++) {
                    int jValue = neighbors[tour[i]][n];
                    int j = reverseTour[jValue];
                    if (j >= i+2) {
                        int nextJ = j+1 < tour.size() ? j+1 : 0;
                        if (tour[i]+tour[i+1] < tour[i] + tour[j]) {
                            break;
                        }
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
        return;
    }
    
    bool kopt2neighbors(vector<int> &tour, vector<int> &reverseTour, const int i) {
        bool didImprove = false;
        for (int n = 1; n < min(200, points); n++) {
            int jValue = neighbors[tour[i]][n];
            int j = reverseTour[jValue];
            if (j >= i+2) {
                int nextJ = j+1 < tour.size() ? j+1 : 0;
                if (dist(tour[i], tour[i+1]) < dist(tour[i], tour[j])) {
                    break;
                }
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
    
    
    vector<int> kopt3(vector<int> &tour) {
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
                            reverse(&tour[i]+1, &tour[j]+1);
                            reverse(&tour[j]+1, &tour[k]+1);
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
    
    void kopt3neighbors(vector<int> &tour) {
        vector<int> reverseTour(tour.size(), 0);
        for (auto i = 0; i < tour.size(); i++) {
            reverseTour[tour[i]] = i;
        }
        int neighborLimit = 20;
        
        bool didImprove = true;
        while (didImprove) {
            didImprove = false;
            for (int i = 0; i < tour.size(); i++) {
                for (int iNeighbor = 1; iNeighbor < min(neighborLimit, (int)tour.size()); iNeighbor++) {
                    int j = reverseTour[neighbors[tour[i]][iNeighbor]];
                    for (int jNeighbor = 1; jNeighbor < min(neighborLimit, (int)tour.size()); jNeighbor++) {
                        int k = reverseTour[neighbors[tour[i+1]][jNeighbor]];
                        if (i < j && j < k) {
                            if (timeIsRunningOut()) {
                                cerr << "Time out" << endl;
                                return;
                            }
                            int nextK = k+1 < tour.size() ? k+1 : 0;
                            double prevDistance2 = dist(tour[i], tour[i+1]) + dist(tour[j], tour[j+1]) + dist(tour[k], tour[nextK]);
                            double afterDistance2 = dist(tour[i], tour[j]) + dist(tour[i+1], tour[k]) + dist(tour[j+1], tour[nextK]);
                            if (afterDistance2 < prevDistance2) {
                                reverse(&tour[i]+1, &tour[j]+1);
                                reverse(&tour[j]+1, &tour[k]+1);
                                vector<int> reverseTour(tour.size(), 0);
                                for (auto x = i+1; x < k+1; x++) {
                                    reverseTour[tour[x]] = x;
                                }
                                didImprove = true;
                            }
                        }
                    }
                }
            }
        }
    }
};


int main(int argc, char **argv) {
#ifdef DEBUG
    auto instance = TravelingSalesmanProblem::createFromFile();
#else
    auto instance = TravelingSalesmanProblem::createFromStdin();
#endif
    auto greedyTour = instance->greedy();
    instance->kopt3(greedyTour);
    vector<int> greedyReverseTour(greedyTour.size());
    for (int i = 0; i < greedyTour.size(); i++) {
        greedyReverseTour[greedyTour[i]] = i;
    }
    vector<int> minimumTour = greedyTour;
    double minimumDistance = instance->tourDistance(greedyTour);
//    while (!timeIsRunningOut()) {
//        instance->kopt2(minimumTour);
//        instance->kopt3(minimumTour);
////                instance->kopt2neighbors(tour);
////        auto tour = greedyTour;
////        auto reverseTour = greedyReverseTour;
////        int iteration = 800;
////        while(!timeIsRunningOut()) {
////            iteration--;
////            if (!instance->kopt2(tour, rand() % instance->points)) {
////                if (iteration <= 0) {
////                    break;
////                }
////            } else {
////                iteration = 800;
////            }
////        }
////        double distance = instance->tourDistance(tour);
////        if (distance < minimumDistance) {
////            minimumDistance = distance;
////            minimumTour = tour;
////        }
//        //        random_shuffle(tour.begin(), tour.end());
//    }
    

    for (auto i: minimumTour) {
        cout << i << endl;
    }
    cerr << "Tour length: " << instance->tourDistance(minimumTour) << endl;
    
    auto endTime = chrono::high_resolution_clock::now();
    auto duration = (endTime - startTime);
    cerr << "Time: " << duration.count() / 1e6 << endl;
    delete instance;
    return 0;
}