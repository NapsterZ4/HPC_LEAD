/**
 * LEAD University
 * Data Science Program
 * BCD-9218: Parallel and Distributed Computing
 * Instructor Diego Jimenez, Eng. (diego.jimenez@ulead.ac.cr)
 *  OpenMP parallel shear sort.
 */

#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <math.h>
#include "timer.h"
#include "io.h"

#define MAX_VALUE 10000

void bubbleSort(int **A, int row, int M){
#pragma omp parallel
    {
#pragma omp parallel for
        for (int i = 0; i < M - 1; ++i) {
            for (int j = 0; j < M - i - 1; ++j) {
                if (A[row][j] > A[row][j + 1]) {
                    int temp = A[row][j];
                    A[row][j] = A[row][j + 1];
                    A[row][j + 1] = temp;
                }
            }
        }
    }
}

void bubbleSortReverse(int **A, int row, int M){
#pragma omp parallel
    {
#pragma omp parallel for
        for (int i = 0; i < M - 1; ++i) {
            for (int j = 0; j < M - i - 1; ++j) {
                if (A[row][j] < A[row][j + 1]) {
                    int temp = A[row][j];
                    A[row][j] = A[row][j + 1];
                    A[row][j + 1] = temp;
                }
            }
        }
    }
}

void bubbleSortCols(int **A, int col, int M){
#pragma omp parallel
    {
#pragma omp parallel for
        for (int i = 0; i < M - 1; ++i) {
            for (int j = 0; j < M - i - 1; ++j) {
                if (A[j][col] > A[j + 1][col]) {
                    int temp = A[j][col];
                    A[j][col] = A[j + 1][col];
                    A[j + 1][col] = temp;
                }
            }
        }
    }
}

// Shear sort function
void shearSort(int **A, int M) {
    int n = (int) ceil(log(M));

//Parallel region
#pragma omp parallel
    {
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < M; ++j) {
                if (j % 2 == 0) { //Sort odd row
                    bubbleSort(A, j, M);
                } else { //Sort even row
                    bubbleSortReverse(A, j, M);
                }
            }
            //Sort columns
            for (int j = 0; j < M; ++j) {
                bubbleSortCols(A, j, M);
            }
        }

        for (int i = 0; i < M; ++i) { //Last phase of algorithm
            if (i % 2 == 0) {
                bubbleSort(A, i, M);
            } else {
                bubbleSortReverse(A, i, M);
            }
        }
    }
}

// Main method
int main(int argc, char *argv[]) {
    int N, M;
    int **A;
    double elapsedTime;

    // checking parameters
    if (argc != 2 && argc != 3) {
        cout << "Parameters: <N> [<file>]" << endl;
        return 1;
    }

    N = atoi(argv[1]);
    M = (int) sqrt(N);

    if (N != M * M) {
        cout << "N has to be a perfect square!" << endl;
        exit(1);
    }

    // allocating matrix A
    A = new int *[M];
    for (int i = 0; i < M; ++i) {
        A[i] = new int[M];
    }

    // reading files (optional)
    if (argc == 3) {
        readMatrixFile(A, M, argv[2]);
    } else {
        srand48(time(NULL));
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                A[i][j] = lrand48() % MAX_VALUE;
            }
        }
    }

    // starting timer
    timerStart();

    // calling shear sort function
    shearSort(A, M);

    // testing the results is correct
    if (argc == 3) {
        printMatrix(A, M);
    }

    // stopping timer
    elapsedTime = timerStop();

    cout << "Duration time: " << elapsedTime << " seconds" << std::endl;

    // releasing memory
    for (int i = 0; i < M; ++i) {
        delete[] A[i];
    }
    delete[] A;

    return 0;
}
