#include <iostream>
#include <chrono>
#include "io.h"

using namespace std;
using namespace std::chrono;

void strassen(int **A, int **B, int **C, int N){


}

int main(int argc, char* argv[]) {
    int N;
    int **A, **B, **C;

    if (argc != 2 && argc != 4){
        cout << " error, Parameters: <Number1> <Number2> <Number3>" << endl;
        return 1;
    }

    N = atoi(argv[1]);

    // Matrix
    A = new int*[N];
    B = new int*[N];
    C = new int*[N];

    for (int i = 0; i < N; ++i) {
        A[i] = new int[N];
        B[i] = new int[N];
        C[i] = new int[N];
    }

    if(argc == 4){
        readMatrixFile(A,N,argv[2]);
        readMatrixFile(B,N,argv[3]);
    }

    // Timestart without function
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // Algoritmo
    strassen(A, B, C, N);

    if (argc == 4){
        printMatrix(C, N);
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    duration <double> totalTime = t2 - t1;
    cout << "Processing time: " << totalTime.count() << " seconds.";

    // Releasing memory
    for (int i=0; i<N; i++) {
        delete [] A[i];
        delete [] B[i];
        delete [] C[i];
    }
    delete [] A;
    delete [] B;
    delete [] C;

    return 0;
}