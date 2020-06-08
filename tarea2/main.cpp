#include <iostream>
#include <chrono>
#include "io.h"
#include "omp.h"

using namespace std;
using namespace std::chrono;

//For 2*2 matrix strassen algorithm
void strassen22(int **A, int **B, int **C) {
    int M1, M2, M3, M4, M5, M6, M7;

    M1 = (A[0][0] + A[1][1]) * (B[0][0] + B[1][1]);
    M2 = (A[1][0] + A[1][1]) * B[0][0];
    M3 = A[0][0] * (B[0][1] - B[1][1]);
    M4 = A[1][1] * (B[1][0] - B[0][0]);
    M5 = (A[0][0] + A[0][1]) * B[1][1];
    M6 = (A[1][0] - A[0][0]) * (B[0][0] + B[0][1]);
    M7 = (A[0][1] - A[1][1]) * (B[1][0] + B[1][1]);

    C[0][0] = M1 + M4 - M5 + M7;
    C[0][1] = M3 + M5;
    C[1][0] = M2 + M4;
    C[1][1] = M1 - M2 + M3 + M6;
}

void strassen(int **A, int **B, int **C, int N){
    //Declare sub-matrix from A and B for multiple position numbers
    int M1[N][N], M2[N][N], M3[N][N], M4[N][N], M5[N][N], M6[N][N], M7[N][N];
    int A11[N][N], A12[N][N], A21[N][N], A22[N][N], B11[N][N], B12[N][N], B21[N][N], B22[N][N];
    int C11[N][N], C12[N][N], C21[N][N], C22[N][N];

    for(int i=0; i<N/2; ++i) {
        for(int j=0; j<N/2; ++j) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j+N/2];
            A21[i][j] = A[i+N/2][j];
            A22[i][j] = A[i+N/2][j+N/2];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j+N/2];
            B21[i][j] = B[i+N/2][j];
            B22[i][j] = B[i+N/2][j+N/2];
        }
    }

    //Calculate M1
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            M1[i][j] = (A11[i][j] + A22[i][j]) * (B11[i][j] + B22[i][j]);
        }
    }

    //Calculate M2
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            M2[i][j] = (A21[i][j] + A22[i][j]) * (B11[i][j]);
        }
    }

    //Calculate M3
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            M3[i][j] = (A11[i][j]) * (B12[i][j] - B22[i][j]);
        }
    }

    //Calculate M4
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            M4[i][j] = (A22[i][j]) * (B21[i][j] - B11[i][j]);
        }
    }

    //Calculate M5
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            M5[i][j] = (A11[i][j] + A12[i][j]) * (B22[i][j]);
        }
    }

    //Calculate M6
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            M6[i][j] = (A21[i][j] - A11[i][j]) * (B11[i][j] + B12[i][j]);
        }
    }

    //Calculate M7
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            M7[i][j] = (A12[i][j] - A22[i][j]) * (B21[i][j] + B22[i][j]);
        }
    }

    //Calculate C1
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            C11[i][j] = M1[i][j] + M4[i][j] + M5[i][j] + M7[i][j];
        }
    }

    //Calculate C2
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            C12[i][j] = M3[i][j] + M5[i][j];
        }
    }

    //Calculate C3
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            C21[i][j] = M2[i][j] + M4[i][j];
        }
    }

    //Calculate C4
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
            C22[i][j] = M1[i][j] + M2[i][j] + M3[i][j] + M6[i][j];
        }
    }

    //Insert in C matrix
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/ 2; ++j) {
            C[i][j] = C11[i][j];
            C[i][j + N/2] = C12[i][j];
            C[i + N/2][j] = C21[i][j];
            C[i + N/2][j + N/2] = C22[i][j];
        }
    }
}

int main(int argc, char* argv[]) {
    int N;
    int **A, **B, **C;

    if (argc != 2 && argc != 4){
        cout << " error, Parameters: <N> <FileA Matrix> <FileB Matrix>" << endl;
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

    // Algorithm
    if (N == 2){
        strassen22(A,B,C);
    } else {
        strassen(A, B, C, N);
    }

    // Impress matrix C
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
