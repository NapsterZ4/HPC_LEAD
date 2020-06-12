#include <iostream>
#include <stdlib.h>
#include <chrono>
#include "io.h"
#include <omp.h>

using namespace std;
using namespace std::chrono;

void strassen(int **A, int **B, int **C, int N, int **A11, int **A12, int **A21, int **A22, int **B11, int **B12, int **B21,
            int **B22, int **M1, int **M2, int **M3, int **M4, int **M5, int **M6, int **M7, int **C11, int **C12,
            int **C21, int **C22, int **AA1, int **AA2, int **AA3, int **AA4, int **AA5, int **BB1, int **BB2, int **BB3,
            int **BB4, int **BB5) {

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
    // Alternatives sub-matrix's
    int **A11, **A12, **A21, **A22, **B11, **B12, **B21, **B22;
    int **M1, **M2, **M3, **M4, **M5, **M6, **M7;
    int **C11, **C12, **C21, **C22;
    int **AA1, **AA2, **AA3, **AA4, **AA5;
    int **BB1, **BB2, **BB3, **BB4, **BB5;

    if (argc != 2 && argc != 4){
        cout << " error, Parameters: <N> <FileA Matrix> <FileB Matrix>" << endl;
        return 1;
    }

    N = atoi(argv[1]);

    // Matrix
    A = new int*[N];
    B = new int*[N];
    C = new int*[N];
    //Sub matrix's
    A11 = new int*[N];
    A12 = new int*[N];
    A21 = new int*[N];
    A22 = new int*[N];
    B11 = new int*[N];
    B12 = new int*[N];
    B21 = new int*[N];
    B22 = new int*[N];
    M1 = new int*[N];
    M2 = new int*[N];
    M3 = new int*[N];
    M4 = new int*[N];
    M5 = new int*[N];
    M6 = new int*[N],
    M7 = new int*[N];
    C11 = new int*[N];
    C12 = new int*[N];
    C21 = new int*[N];
    C22 = new int*[N];
    AA1 = new int*[N];
    AA2 = new int*[N];
    AA3 = new int*[N];
    AA4 = new int*[N];
    AA5 = new int*[N];
    BB1 = new int*[N];
    BB2 = new int*[N];
    BB3 = new int*[N];
    BB4 = new int*[N];
    BB5 = new int*[N];

    for (int i = 0; i < N; ++i) {
        // Matrix
        A[i] = new int[N];
        B[i] = new int[N];
        C[i] = new int[N];
        //Sub matrix's
        A11[i] = new int[N];
        A12[i] = new int[N];
        A21[i] = new int[N];
        A22[i] = new int[N];
        B11[i] = new int[N];
        B12[i] = new int[N];
        B21[i] = new int[N];
        B22[i] = new int[N];
        M1[i] = new int[N];
        M2[i] = new int[N];
        M3[i] = new int[N];
        M4[i] = new int[N];
        M5[i] = new int[N];
        M6[i] = new int[N],
        M7[i] = new int[N];
        C11[i] = new int[N];
        C12[i] = new int[N];
        C21[i] = new int[N];
        C22[i] = new int[N];
        AA1[i] = new int[N];
        AA2[i] = new int[N];
        AA3[i] = new int[N];
        AA4[i] = new int[N];
        AA5[i] = new int[N];
        BB1[i] = new int[N];
        BB2[i] = new int[N];
        BB3[i] = new int[N];
        BB4[i] = new int[N];
        BB5[i] = new int[N];
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
    strassen(A, B, C, N, A11, A12, A21, A22, B11, B12, B21, B22, M1, M2, M3, M4, M5, M6, M7, C11, C12, C21, C22,
             AA1, AA2, AA3, AA4, AA5, BB1, BB2, BB3, BB4, BB5);

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
        delete [] A11[i];
        delete [] A12[i];
        delete [] A21[i];
        delete [] A22[i];
        delete [] B11[i];
        delete [] B12[i];
        delete [] B21[i];
        delete [] B22[i];
        delete [] M1[i];
        delete [] M2[i];
        delete [] M3[i];
        delete [] M4[i];
        delete [] M5[i];
        delete [] M6[i];
        delete [] M7[i];
        delete [] C11[i];
        delete [] C12[i];
        delete [] C21[i];
        delete [] C22[i];
        delete [] AA1[i];
        delete [] AA2[i];
        delete [] AA3[i];
        delete [] AA4[i];
        delete [] AA5[i];
        delete [] BB1[i];
        delete [] BB2[i];
        delete [] BB3[i];
        delete [] BB4[i];
        delete [] BB5[i];
    }
    delete [] A;
    delete [] B;
    delete [] C;
    delete [] A11;
    delete [] A12;
    delete [] A21;
    delete [] A22;
    delete [] B11;
    delete [] B12;
    delete [] B21;
    delete [] B22;
    delete [] M1;
    delete [] M2;
    delete [] M3;
    delete [] M4;
    delete [] M5;
    delete [] M6;
    delete [] M7;
    delete [] C11;
    delete [] C12;
    delete [] C21;
    delete [] C22;
    delete [] AA1;
    delete [] AA2;
    delete [] AA3;
    delete [] AA4;
    delete [] AA5;
    delete [] BB1;
    delete [] BB2;
    delete [] BB3;
    delete [] BB4;
    delete [] BB5;

    return 0;
}