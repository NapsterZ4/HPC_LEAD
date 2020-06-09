#include <iostream>
#include <stdlib.h>
#include <chrono>
#include "io.h"
#include <omp.h>

using namespace std;
using namespace std::chrono;

void strassen(int **A, int **B, int **C, int N){
    //Declare sub-matrix from A and B for multiple position numbers and alternatives matrix
    int M1[N][N], M2[N][N], M3[N][N], M4[N][N], M5[N][N], M6[N][N], M7[N][N];
    int A11[N][N], A12[N][N], A21[N][N], A22[N][N], B11[N][N], B12[N][N], B21[N][N], B22[N][N];
    int C11[N][N], C12[N][N], C21[N][N], C22[N][N];
    int AA1[N][N], AA2[N][N], AA3[N][N], AA4[N][N], AA5[N][N];
    int BB1[N][N], BB2[N][N], BB3[N][N], BB4[N][N], BB5[N][N];

    for (int i = 0; i < N / 2; ++i) {
        for (int j = 0; j < N / 2; ++j) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + N / 2];
            A21[i][j] = A[i + N / 2][j];
            A22[i][j] = A[i + N / 2][j + N / 2];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + N / 2];
            B21[i][j] = B[i + N / 2][j];
            B22[i][j] = B[i + N / 2][j + N / 2];
        }
    }

    //First parallel region
#pragma omp parallel
    {
#pragma omp task
        //Calculate M1
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                AA1[i][j] = (A11[i][j] + A22[i][j]); //add
            }
        }

#pragma omp task
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                BB1[i][j] = (B11[i][j] + B22[i][j]); //add
            }
        }

#pragma omp task
        //Insert M1 multiply
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                M1[i][j] = 0;
                for (int k = 0; k < N / 2; ++k) {
                    M1[i][j] = M1[i][j] + AA1[i][k] * BB1[k][j];
                }
            }
        }

#pragma omp task
        //Calculate M2
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                AA2[i][j] = (A21[i][j] + A22[i][j]); //add
            }
        }

#pragma omp task
        //Insert M2 multiply
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                M2[i][j] = 0;
                for (int k = 0; k < N / 2; ++k) {
                    M2[i][j] = M2[i][j] + AA2[i][k] * B11[k][j];
                }
            }
        }

#pragma omp task
        //Calculate M3
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                BB2[i][j] = B12[i][j] - B22[i][j]; //sub
            }
        }

#pragma omp task
        //Insert M3 multiply
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                M3[i][j] = 0;
                for (int k = 0; k < N / 2; ++k) {
                    M3[i][j] = M3[i][j] + A11[i][k] * (BB2[k][j]);
                }
            }
        }

#pragma omp task
        //Calculate M4
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                BB3[i][j] = B21[i][j] - B11[i][j]; //sub
            }
        }

#pragma omp task
        //Insert M4 multiply
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                M4[i][j] = 0;
                for (int k = 0; k < N / 2; ++k) {
                    M4[i][j] = M4[i][j] + A22[i][k] * BB3[k][j];
                }
            }
        }

#pragma omp task
        //Calculate M5
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                AA3[i][j] = A11[i][j] + A12[i][j]; //add
            }
        }

#pragma omp task
        //Insert M5 multiply
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                M5[i][j] = 0;
                for (int k = 0; k < N / 2; ++k) {
                    M5[i][j] = M5[i][j] + AA3[i][k] * B22[k][j];
                }
            }
        }

#pragma mp task
        //Calculate M6
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                AA4[i][j] = A21[i][j] - A11[i][j]; //sub
            }
        }

#pragma omp task
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                BB4[i][j] = B11[i][j] + B12[i][j]; //add
            }
        }

#pragma omp task
        //Insert M6 multiply
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                M6[i][j] = 0;
                for (int k = 0; k < N / 2; ++k) {
                    M6[i][j] = M6[i][j] + AA4[i][k] * BB4[k][j];
                }
            }
        }

#pragma omp task
        //Calculate M7
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                AA5[i][j] = A12[i][j] - A22[i][j]; //sub
            }
        }
    }

    //Second parallel region
#pragma omp parallel
    {
#pragma omp task
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                BB5[i][j] = B21[i][j] + B22[i][j]; //add
            }
        }

#pragma omp task
        //Insert M6 multiply
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                M7[i][j] = 0;
                for (int k = 0; k < N / 2; ++k) {
                    M7[i][j] = M7[i][j] + AA5[i][k] * BB5[k][j];
                }
            }
        }

#pragma omp task
        //Calculate C1
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                C11[i][j] = M1[i][j] + M4[i][j] - M5[i][j] + M7[i][j];
            }
        }

#pragma omp task
        //Calculate C2
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                C12[i][j] = M3[i][j] + M5[i][j];
            }
        }

#pragma omp task
        //Calculate C3
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                C21[i][j] = M2[i][j] + M4[i][j];
            }
        }

#pragma omp task
        //Calculate C4
        for (int i = 0; i < N / 2; ++i) {
            for (int j = 0; j < N / 2; ++j) {
                C22[i][j] = M1[i][j] - M2[i][j] + M3[i][j] + M6[i][j];
            }
        }
    }

#pragma omp taskwait
    //Insert in C matrix
    for (int i = 0; i < N/2; ++i) {
        for (int j = 0; j < N/2; ++j) {
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
    } else {
        //Matrix random generator
        matrixRandomComplete(A,B,N);
    }

    // Timestart without function
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // Algorithm
    strassen(A, B, C, N);

    // Impress matrix C
    if (argc == 4){
        printMatrix(C, N);
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    duration <double> totalTime = t2 - t1;
    cout << "Processing time: " << totalTime.count() << " seconds.\n";

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