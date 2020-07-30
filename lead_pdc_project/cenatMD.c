/**
 * Costa Rica National High Technology Center
 * Advanced Computing Laboratory
 * MPI molecular dynamics 
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// General constants
#define TAG 7
#define CONSTANT 777

// Particle interaction constants
#define A 2.0
#define B 1.0
#define MASS 4
#define DELTA 0.1
#define MAX_DISTANCE 4.0
#define LIMIT 60.0
#define INIT_VELOCITY 0.0

// Random initialization constants
#define POSITION 0
#define POSITION_CENTER 1
#define POSITION_X 2
#define VELOCITY 3

// Structure for shared properties of a particle
struct particle {
    double x;
    double y;
    double z;
    double mass;
    double fx;
    double fy;
    double fz;
    double aX;
    double aY;
    double aZ;
    double vX;
    double vY;
    double vZ;
};

// Heading for auxiliar functions
double randomValue(int type, int rank, int p);

void printFile(struct particle *shared, int iteration, int rank, int limit);

void printInitFile(struct particle *shared, int limit);

void interact(struct particle *source, struct particle *destiny);

void evolve(struct particle *source, struct particle *destiny, int limitSource, int limitDestiny);

void merge(struct particle *first, struct particle *second, int limit);

void updateProperties(struct particle *shared, int limit);

void cleanForces(struct particle *shared, int limit);

double randomRange(double min, double max);

// Main function
int main(int argc, char **argv) {
    int myRank;                                    // Rank of process
    int p;                                        // Number of processes
    int n;                                        // Base number of particles
    int iterations;                                    // Number of iterations
    int previous;                                    // Previous process in the ring
    int next;                                    // Next process in the ring
    int tag = TAG;                                    // Tag for message
    int number, maxNumber, foreignNumber;                        // Number of particles
    struct particle *locals;                            // Array of local particles
    struct particle *foreigners;                            // Array of foreign particles
    MPI_Status status;                                // Return status for receive
    int i, j, rounds, origin, sender, dump_flag, init_flag;
    double start_time, end_time;

    // checking the number of parameters
    if (argc < 5) {
        printf("Error, usage %s <particles per rank> <iterations> <dump flag> <init flag>\n", argv[0]);
        exit(1);
    }

    // getting base number of particles per rank
    n = atoi(argv[1]);
    iterations = atoi(argv[2]);
    dump_flag = atoi(argv[3]);
    init_flag = atoi(argv[4]);

    // initializing MPI structures and checking for non-parity of p
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if (p % 2 == 0) {
        p = p - 1;
        if (myRank == p) {
            MPI_Finalize();
            return 0;
        }
    }

    // setting random seed
    srand48(myRank + myRank * CONSTANT);

    //TO DO: compute previous and next processes
    next = (myRank + 1) % p;
    previous = (myRank + p - 1) % p;

    number = n;
    maxNumber = n;

    // acquiring memory for particle arrays
    locals = malloc(number * sizeof(struct particle));
    foreigners = malloc(maxNumber * sizeof(struct particle));

    // initializing particle positions
    if (init_flag) {
        for (j = 0; j < number; j++) {
            locals[j].x = (double) ((myRank + j) % (int) LIMIT);
            locals[j].y = (double) ((myRank + 2 * j) % (int) LIMIT);
            locals[j].z = (double) ((number + j) % (int) LIMIT);
            locals[j].fx = 0.0;
            locals[j].fy = 0.0;
            locals[j].fz = 0.0;
            locals[j].mass = MASS;
            locals[j].vX = INIT_VELOCITY;
            locals[j].vY = INIT_VELOCITY;
            locals[j].vZ = INIT_VELOCITY;
        }
    } else {
        // random initialization of local particle array
        for (j = 0; j < number; j++) {
            if (j % 2 == 0) {
                locals[j].x = randomValue(POSITION_X, myRank, p);
                locals[j].y = randomValue(POSITION, myRank, p);
                locals[j].z = randomValue(POSITION, myRank, p);
            } else {
                locals[j].x = randomValue(POSITION_X, myRank, p);
                locals[j].y = randomValue(POSITION_CENTER, myRank, p);
                locals[j].z = randomValue(POSITION_CENTER, myRank, p);
            }
            locals[j].fx = 0.0;
            locals[j].fy = 0.0;
            locals[j].fz = 0.0;
            locals[j].mass = MASS;
            locals[j].vX = randomValue(VELOCITY, myRank, p);
            locals[j].vY = randomValue(VELOCITY, myRank, p);
            locals[j].vZ = randomValue(VELOCITY, myRank, p);
        }
    }


    start_time = MPI_Wtime();

    // executing iterations
    for (i = 1; i <= iterations; i++) {

        // cleaning forces in the particles
        cleanForces(locals, number);
        foreignNumber = n;

        //TO DO: sending the local particles to the next processor, receiving the incoming foreign particle set and update both of them
        MPI_Send(locals, number * (sizeof(struct particle)) / sizeof(double), MPI_DOUBLE,
                 next, tag, MPI_COMM_WORLD);

        MPI_Recv(foreigners, number * (sizeof(struct particle)) / sizeof(double), MPI_DOUBLE, previous, tag,
                 MPI_COMM_WORLD, &status);


        evolve(locals, foreigners, number, foreignNumber);

        //TO DO: running the algorithm for (p-1)/2 rounds. REMEMBER: call evolve function
        int ring_iterations = (p - 1) / 2;
        for (int i = 0; i < ring_iterations; ++i) {
            MPI_Send(locals, number * (sizeof(struct particle)) / sizeof(double), MPI_DOUBLE,
                     next, tag, MPI_COMM_WORLD);

            MPI_Recv(foreigners, number * (sizeof(struct particle)) / sizeof(double), MPI_DOUBLE, previous, tag,
                     MPI_COMM_WORLD, &status);

            evolve(locals, foreigners, number, foreignNumber);
        }

        //TO DO: sending the particles to the origin
//        origin = (myRank + 1) % ((p - 1) / 2);
//        origin = (myRank + (p - 1) / 2) % p;
        origin = (myRank + 1 + ring_iterations) % p;
        int last = (myRank + ring_iterations) % p;

        MPI_Send(foreigners, number * (sizeof(struct particle)) / sizeof(double), MPI_DOUBLE,
                 origin, tag, MPI_COMM_WORLD);

        MPI_Recv(locals, number * (sizeof(struct particle)) / sizeof(double), MPI_DOUBLE,
                 last, tag, MPI_COMM_WORLD, &status);

        //TO DO: receiving the incoming particles and merging them with the local set, interacting the local set

        merge(locals, foreigners, number);
        evolve(locals, locals, number, number);

        // computing new velocity, acceleration and position
        updateProperties(locals, number);

        // dumping particle positions on a file
        if (dump_flag && (i % 100 == 0)) {
            if (myRank == 0) {
                int k;
                printFile(locals, i, myRank, number);
                for (k = 1; k < p; k++) {
                    MPI_Recv(foreigners, number * (sizeof(struct particle)) / sizeof(double), MPI_DOUBLE, k, tag,
                             MPI_COMM_WORLD, &status);
                    printFile(foreigners, i, myRank, number);
                }
            } else {
                MPI_Send(locals, number * (sizeof(struct particle)) / sizeof(double), MPI_DOUBLE, 0, tag,
                         MPI_COMM_WORLD);
            }
        }
    }

    end_time = MPI_Wtime();
    if (myRank == 0)
        printf("Elapsed Time: %f seconds\n", end_time - start_time);

    // disposing memory for particle arrays
    free(locals);
    free(foreigners);

    // finalizing MPI structures
    MPI_Finalize();

    return 0;
}

// Function for random value generation
double randomValue(int type, int rank, int p) {
    double value;
    switch (type) {
        case POSITION:
            value = (double) drand48() * (double) LIMIT;
            break;
        case POSITION_X:
            value = randomRange((double) rank * (double) LIMIT / (double) p,
                                (double) (rank + 1) * (double) LIMIT / (double) p);
            break;
        case POSITION_CENTER:
            value = randomRange((double) LIMIT / 3.0, (double) LIMIT * 2.0 / 3.0);
            break;
        case VELOCITY:
            value = INIT_VELOCITY;
            break;
        default:
            value = 1.1;
    }
    return value;
}

// Auxiliary function to generate a random number in a given range
double randomRange(double min, double max) {
    return min + (double) drand48() * (max - min);
}

// Auxiliary function to print particle positions on a file
void printFile(struct particle *shared, int iteration, int rank, int limit) {
    FILE *output;
    int j;
    char filename[64];
    sprintf(filename, "cenatMD_%d.txt", iteration);


    output = fopen(filename, "a+r");

    int c = fgetc(output);
    if (c == EOF) {
        fprintf(output, "x,y,z\n");
    } else {
        ungetc(c, output);
    }


    for (j = 0; j < limit; j++) {
        fprintf(output, "%7.12f,%7.12f,%7.12f\n", shared[j].x, shared[j].y, shared[j].z);

    }

    fclose(output);

}


void printInitFile(struct particle *shared, int limit) {
    FILE *output;
    int j;
    char filename[64];
    sprintf(filename, "cenatMD_init.txt");


    output = fopen(filename, "a+r");

    int c = fgetc(output);
    if (c == EOF) {
        fprintf(output, "x,y,z\n");
    } else {
        ungetc(c, output);
    }


    for (j = 0; j < limit; j++) {
        fprintf(output, "%7.12f,%7.12f,%7.12f\n", shared[j].x, shared[j].y, shared[j].z);

    }

    fclose(output);

}


// Function for computing interaction among two particles
// There is an extra test for interaction of identical particles, in which case there is no effect over the destiny
void interact(struct particle *first, struct particle *second) {
    double rx, ry, rz, r, fx, fy, fz, f;

    // computing base values
    rx = first->x - second->x;
    ry = first->y - second->y;
    rz = first->z - second->z;
    r = sqrt(rx * rx + ry * ry + rz * rz);

    if (r == 0.0 || r >= MAX_DISTANCE)
        return;

    f = A / pow(r, 12) - B / pow(r, 6);
    fx = f * rx / r;
    fy = f * ry / r;
    fz = f * rz / r;

    //printf("Datos de fuerza: %.12f, %.12f, %.12f\n",fx,fy,fz);

    // updating destiny's structure
    second->fx = second->fx + fx;
    second->fy = second->fy + fy;
    second->fz = second->fz + fz;

    // updating sources's structure
    first->fx = first->fx - fx;
    first->fy = first->fy - fy;
    first->fz = first->fz - fz;
}

// Function for computing interaction among two particle arrays
void evolve(struct particle *first, struct particle *second, int limitFirst, int limitSecond) {
    int j, k;

    for (j = 0; j < limitFirst; j++) {
        for (k = j + 1; k < limitSecond; k++) {
            interact(&first[j], &second[k]);
        }
    }
}

// Function to merge two particle arrays
// Permanent changes reside only in first array
void merge(struct particle *first, struct particle *second, int limit) {
    int j;

    for (j = 0; j < limit; j++) {story |
        first[j].fx += second[j].fx;
        first[j].fy += second[j].fy;
        first[j].fz += second[j].fz;
    }
}

// Function for updating velocity, acceleration and new position
void updateProperties(struct particle *shared, int limit) {
    int j;
    double ax, ay, vx, vy, x, y;

    for (j = 0; j < limit; j++) {
        shared[j].aX = shared[j].fx / MASS;
        shared[j].aY = shared[j].fy / MASS;
        shared[j].aZ = shared[j].fz / MASS;
        shared[j].vX = shared[j].vX + shared[j].aX * DELTA;
        shared[j].vY = shared[j].vY + shared[j].aY * DELTA;
        shared[j].vZ = shared[j].vZ + shared[j].aZ * DELTA;
        shared[j].x = shared[j].x + shared[j].vX * DELTA;
        shared[j].y = shared[j].y + shared[j].vY * DELTA;
        shared[j].z = shared[j].z + shared[j].vZ * DELTA;

        // reflecting particules
        while (shared[j].x > LIMIT || shared[j].x < 0.0) {
            if (shared[j].x > LIMIT)
                shared[j].x = shared[j].x - (shared[j].x - LIMIT);
            else if (shared[j].x < 0.0)
                shared[j].x *= -1;
            shared[j].vX *= -1;
        }
        while (shared[j].y > LIMIT || shared[j].y < 0.0) {
            if (shared[j].y > LIMIT)
                shared[j].y = shared[j].y - (shared[j].y - LIMIT);
            else if (shared[j].y < 0.0)
                shared[j].y *= -1;
            shared[j].vY *= -1;
        }
        while (shared[j].z > LIMIT || shared[j].z < 0.0) {
            if (shared[j].z > LIMIT)
                shared[j].z = shared[j].z - (shared[j].z - LIMIT);
            else if (shared[j].z < 0.0)
                shared[j].z *= -1;
            shared[j].vZ *= -1;
        }
    }
}

// Function to clean forces in particles
void cleanForces(struct particle *shared, int limit) {
    int j;

    for (j = 0; j < limit; j++) {
        shared[j].fx = 0.0;
        shared[j].fy = 0.0;
        shared[j].fz = 0.0;
    }
}