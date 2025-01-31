#include <mpi.h>
#include <stdio.h>
#include <string.h>

void BCS221088_Input(int *, int *, double **, double **, double **, int);
void BCS221088_Multiplication(int, int, int, int, double *, double *, double *);
void BCS221088_Print(int, int, int, double *);
void BCS221088_Mytime(int, double);


int main() {
    int my_rank, comm_size;
    int matrixRows = 0, matrixCols = 0;
    double *fullMatrix = NULL, *vector = NULL;
    double *localMatrixSegment = NULL, *localVectorSegment = NULL;
    double startTimestamp, endTimestamp;

    // Initialize MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // Synchronize processes and record start time
    MPI_Barrier(MPI_COMM_WORLD);
    startTimestamp = MPI_Wtime();

    // Read input and distribute data among processes
    BCS221088_Input(&matrixRows, &matrixCols, &fullMatrix, &vector, &localMatrixSegment, my_rank);

    // Allocate memory for the local result
    int segmentRows = matrixRows / comm_size;
    localVectorSegment = (double *)malloc(segmentRows * sizeof(double));

    // Perform matrix-vector multiplication
    BCS221088_Multiplication(my_rank, comm_size, segmentRows, matrixCols, localMatrixSegment, vector, localVectorSegment);

    // Display the final result
    BCS221088_Print(my_rank, comm_size, matrixRows, localVectorSegment);

    // Record end time and calculate elapsed time
    endTimestamp = MPI_Wtime();
    elapsed_local = endTimestamp - startTimestamp;
    BCS221088_Mytime(my_rank, elapsed_local);

    // Free allocated memory
    if (my_rank == 0) {
        free(fullMatrix);
        free(vector);
    }
    free(localMatrixSegment);
    free(localVectorSegment);

    // Finalize MPI
    MPI_Finalize();
    return 0;
}

// Function to read input and distribute data among processes
void BCS221088_Input(int *rowCount, int *colCount, double **matrix, double **vector, double **localMatrix, int my_rank) {
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if (my_rank == 0) {
        printf("Enter the number of rows: ");
        scanf("%d", rowCount);

        printf("Enter the number of columns: ");
        scanf("%d", colCount);

        *vector = (double *)malloc((*colCount) * sizeof(double));
        *matrix = (double *)malloc((*rowCount) * (*colCount) * sizeof(double));

        printf("Enter matrix values:\n");
        for (int rIndex = 0; rIndex < *rowCount; rIndex++) {
            for (int cIndex = 0; cIndex < *colCount; cIndex++) {
                printf("Matrix[%d][%d]: ", rIndex, cIndex);
                scanf("%lf", (*matrix + rIndex * (*colCount) + cIndex));
            }
        }

        printf("Enter vector values:\n");
        for (int vIndex = 0; vIndex < *colCount; vIndex++) {
            printf("Vector[%d]: ", vIndex);
            scanf("%lf", (*vector + vIndex));
        }
    }
    if (my_rank != 0){
        *vector = (double *)malloc((*colCount) * sizeof(double));
    }
    int segmentRows = (*rowCount) / comm_size;
    *localMatrix = (double *)malloc(segmentRows * (*colCount) * sizeof(double));
    
    MPI_Bcast(rowCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(colCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(*vector, *colCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(*matrix, segmentRows * (*colCount), MPI_DOUBLE, *localMatrix, segmentRows * (*colCount), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

// Function to perform matrix-vector multiplication
void BCS221088_Multiplication(int my_rank, int comm_size, int rowsPerProcess, int colCount, double *localMatrix, double *vector, double *localResult) {
    for (int row = 0; row < rowsPerProcess; row++) {
        localResult[row] = 0.0;
        for (int col = 0; col < colCount; col++) {
            localResult[row] += localMatrix[row * colCount + col] * vector[col];
        }
    }
}

// Function to display the result
void BCS221088_Print(int my_rank, int comm_size, int totalRows, double *localResult) {
    double *globalResult = NULL;
    int segmentRows = totalRows / comm_size;

    if (my_rank == 0) {
        globalResult = (double *)malloc(totalRows * sizeof(double));
    }

    MPI_Gather(localResult, segmentRows, MPI_DOUBLE, globalResult, segmentRows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("Resultant vector:\n");
        for (int index = 0; index < totalRows; index++) {
            printf("y[%d]: %.2lf\n", index, globalResult[index]);
        }
        free(globalResult);
    }
}

// Function to calculate and display elapsed time
// Function to calculate and display the maximum elapsed time
void BCS221088_Mytime(int my_rank, double elapsed) {
    double maxElapsedTime;

    // Find the maximum elapsed time among all processes
    MPI_Reduce(&elapsedTime, &maxElapsedTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("Maximum time taken by any process: %.6lf seconds\n", maxElapsedTime);
    }
}

