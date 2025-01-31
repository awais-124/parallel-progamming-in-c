#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Function prototypes
void bcs221098_Input(int rank, int comm_size, double** local_A, double* local_x, int* local_m, int* n, MPI_Comm comm);
void bcs221098_Multiplication(double local_A[], double local_x[], double local_y[], int local_m, int n);
void bcs221098_Print(double local_y[], int local_m, int n, MPI_Comm comm, int rank);
double bcs221098_Mytime();

int main(int argc, char** argv) {
    int rank, comm_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int n;           // Number of columns in matrix (and size of vector x)
    int m;           // Number of rows in matrix
    int local_m;     // Number of rows per process

    double* local_A;
    double* local_x;
    double* local_y;

    // Input the matrix and vector dimensions, and allocate memory
    bcs221098_Input(rank, comm_size, &local_A, &local_x, &local_m, &n, MPI_COMM_WORLD);
    local_y = (double*)malloc(local_m * sizeof(double));

    // Start timing
    double start_time = bcs221098_Mytime();

    // Perform matrix-vector multiplication
    bcs221098_Multiplication(local_A, local_x, local_y, local_m, n);

    // Gather and print the result
    bcs221098_Print(local_y, local_m, n, MPI_COMM_WORLD, rank);

    // End timing
    if (rank == 0) {
        double end_time = bcs221098_Mytime();
        printf("Time taken: %f seconds\n", end_time - start_time);
    }

    free(local_A);
    free(local_x);
    free(local_y);

    MPI_Finalize();
    return 0;
}

void bcs221098_Input(int rank, int comm_size, double** local_A, double* local_x, int* local_m, int* n, MPI_Comm comm) {
    if (rank == 0) {
        // Master process inputs the matrix dimensions
        int m;
        printf("Enter the number of rows in the matrix (m): ");
        scanf("%d", &m);

        printf("Enter the number of columns in the matrix (n): ");
        scanf("%d", n);

        *local_m = m / comm_size;

        // Allocate memory for the full matrix and vector
        double* A = (double*)malloc(m * (*n) * sizeof(double));
        double* x = (double*)malloc((*n) * sizeof(double));

        // Input the matrix and vector
        printf("Enter the matrix A (%d x %d):\n", m, *n);
        for (int i = 0; i < m * (*n); i++) {
            scanf("%lf", &A[i]);
        }

        printf("Enter the vector x (%d x 1):\n", *n);
        for (int i = 0; i < *n; i++) {
            scanf("%lf", &x[i]);
        }

        // Scatter matrix and broadcast vector
        local_A = (double)malloc((*local_m) * (*n) * sizeof(double));
        MPI_Scatter(A, (*local_m) * (*n), MPI_DOUBLE, *local_A, (*local_m) * (*n), MPI_DOUBLE, 0, comm);
        local_x = (double)malloc((*n) * sizeof(double));
        MPI_Bcast(x, *n, MPI_DOUBLE, 0, comm);

        free(A);
        free(x);
    } else {
        // Receive the dimensions from the master process
        MPI_Bcast(n, 1, MPI_INT, 0, comm);
        int m;
        MPI_Bcast(&m, 1, MPI_INT, 0, comm);
        *local_m = m / comm_size;

        // Allocate memory for local matrix and vector
        local_A = (double)malloc((*local_m) * (*n) * sizeof(double));
        local_x = (double)malloc((*n) * sizeof(double));

        // Receive scattered matrix and broadcast vector
        MPI_Scatter(NULL, (*local_m) * (*n), MPI_DOUBLE, *local_A, (*local_m) * (*n), MPI_DOUBLE, 0, comm);
        MPI_Bcast(*local_x, *n, MPI_DOUBLE, 0, comm);
    }
}

void bcs221098_Multiplication(double local_A[], double local_x[], double local_y[], int local_m, int n) {
    for (int i = 0; i < local_m; i++) {
        local_y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            local_y[i] += local_A[i * n + j] * local_x[j];
        }
    }
}

void bcs221098_Print(double local_y[], int local_m, int n, MPI_Comm comm, int rank) {
    if (rank == 0) {
        double* y = (double*)malloc(n * sizeof(double));
        MPI_Gather(local_y, local_m, MPI_DOUBLE, y, local_m, MPI_DOUBLE, 0, comm);

        printf("Resulting vector y:\n");
        for (int i = 0; i < n; i++) {
            printf("%f\n", y[i]);
        }

        free(y);
    } else {
        MPI_Gather(local_y, local_m, MPI_DOUBLE, NULL, local_m, MPI_DOUBLE, 0, comm);
    }
}

double bcs221098_Mytime() {
    return MPI_Wtime();
}