/**
 * Costa Rica National High Technology Center
 * Advanced Computing Laboratory
 * Instructor: Diego Jimenez Vargas, Eng.
 * MPI parallel sum program.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define FLAG 7

// Main routine
int main(int argc, char *argv[]) {
    int rank, size, number, counter;

    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // computing sum of all numbers in the system
    srand(time(NULL) * rank);
    number = rand() % size;

    // verifying result
    MPI_Reduce(&number, &counter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("[%d] Total sum: %d\n", rank, counter);
    }

    // parallel sum algorithm
    for (int i = 1; i <= (int) log2(size); ++i) {
        if (rank % (int) pow(2, i) == pow(2, i) - 1) {
            // TO DO: receive number from leaf on spanning tree
            MPI_Recv(&counter, 1, MPI_INT, rank - (int) pow(2, i - 1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            number += counter;
        }

        if (rank % (int) pow(2, i) == pow(2, i) - pow(2, i - 1) - 1) {
            // TO DO: send number to parent on spanning tree
            MPI_Send(&number, 1, MPI_INT,  rank + (int) pow(2, i - 1), 0, MPI_COMM_WORLD);
        }
    }

    if (rank == size - 1) {
        printf("[%d] Total sum: %d\n", rank, number);
    }

    // finalize MPI
    MPI_Finalize();
    return 0;
}
