/**
 * Costa Rica National High Technology Center
 * Advanced Computing Laboratory
 * Instructor:  Diego Jimenez Vargas, Eng.
 * MPI 2D stencil computation.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define FLAG 7
#define HIGH 200
#define LOW 0
#define TEST 1

// Main routine
int main(int argc, char *argv[]) {
    int rank, size, M, N, rows;
    double **data, **next, **tmp, left, right, start_time, end_time;

    // reading command-line arguments
    if (argc < 3) {
        printf("Error, usage %s <iterations> <grid size>\n", argv[0]);
        exit(1);
    }
    M = atoi(argv[1]);
    N = atoi(argv[2]);

    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // allocating memory for data
    rows = N / size;
    data = (double **) malloc((rows + 2) * sizeof(double *));
    next = (double **) malloc((rows + 2) * sizeof(double *));
    for (int i = 0; i < rows + 2; ++i) {
        data[i] = (double *) malloc(N * sizeof(double));
        next[i] = (double *) malloc(N * sizeof(double));
    }

    for (int i = 0; i < rows + 2; ++i) {
        for (int j = 0; j < N; ++j) {
            data[i][j] = LOW;
        }
    }

    // setting hot plates on both extremes
    if (rank == 0) {
        for (int j = 0; j < N; ++j) {
            data[0][j] = HIGH;
        }
    }

    if (rank == (size - 1)) {
        for (int j = 0; j < N; ++j) {
            data[rows + 1][j] = HIGH;
        }
    }

    start_time = MPI_Wtime();

    tmp = 0;
    // stencil loop
    for (int i = 0; i < M; ++i) {
        // TO DO: send ghost cells to neighbors
        // TO DO: receive ghost cells from neighbors

        if (rank == 0) {
            printf("Starting iteration %d\n", i);
            MPI_Send(&data, sizeof(data), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&tmp, sizeof(data), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else if (rank != 0 && rank != (size - 1)){
            MPI_Send(&data, sizeof(data), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&tmp, sizeof(data), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else if (rank == (size - 1)){
            MPI_Send(&data, sizeof(data), MPI_DOUBLE, size - 2, 0, MPI_COMM_WORLD);
            MPI_Recv(&tmp, sizeof(data), MPI_DOUBLE, size - 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // modify internal elements
        for (int i = 1; i < rows + 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                next[i][j] = (data[i][j] + data[i][j - 1] + data[i][j + 1] + data[i - 1][j] + data[i + 1][j]) / 5.0;
            }
        }
        // switch arrays
        tmp = data;
        data = next;
        next = tmp;
    }

    end_time = MPI_Wtime();
    if (rank == 0) {
        printf("Elapsed Time: %f seconds\n", end_time - start_time);
    }

#if TEST
    if(rank == 0){
        for(int j=0; j<N; j++)
            printf("%.2f ",data[1][j]);
        printf("\n");
    }
#endif

    // finalize MPI
    MPI_Finalize();
    return 0;
}