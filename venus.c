#include <mpi.h>
#include <stdio.h>

// Function prototypes
void BCS221050_Input(int my_rank, double **matrix, double **vector, double **local_matrix, int *row_count, int *col_count);
void BCS221050_Multiplication(int col_count, int rows_per_process, double *local_matrix, double *vector, double *local_result);
void BCS221050_Print(int my_rank, double *local_result, int row_count, int comm_size);
void BCS221050_Mytime(int my_rank, double start_time, double end_time);

int main() {
    // Variables for matrix and vector data
    int row_count = 0, col_count = 0;
    double *matrix = NULL, *vector = NULL, *local_matrix = NULL, *local_result = NULL;

    // Variables for timing
    double start_time, end_time;

    // Initialize MPI environment
    // Variables for MPI rank and size
    int my_rank, comm_size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // Synchronize processes and record the start time
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // Input data and distribute among processes
    BCS221050_Input(my_rank, &matrix, &vector, &local_matrix, &row_count, &col_count);

    // Determine rows to process per rank
    int rows_per_process = row_count / comm_size;
    local_result = (double *)malloc(rows_per_process * sizeof(double));
    BCS221050_Multiplication(col_count, rows_per_process, local_matrix, vector, local_result);
    BCS221050_Print(my_rank, local_result, row_count, comm_size);
    end_time = MPI_Wtime();
    BCS221050_Mytime(my_rank, start_time, end_time);

    // Free allocated memory
    if (my_rank == 0) {
        free(matrix);
        free(vector);
    }
    free(local_matrix);
    free(local_result);

    MPI_Finalize();
    return 0;
}

/*
 * Function to handle input and data distribution among processes.
 */
void BCS221050_Input(int my_rank, double **matrix, double **vector, double **local_matrix, int *row_count, int *col_count) {
    
    if (my_rank == 0) {
        printf("Please enter the number of rows and columns of the matrix (space-separated): ");
        scanf("%d %d", row_count, col_count);

        *vector = (double *)malloc((*col_count) * sizeof(double));
        *matrix = (double *)malloc((*row_count) * (*col_count) * sizeof(double));

        printf("Enter the vector values (space-separated):\n");
        for (int i = 0; i < *col_count; i++) {
            scanf("%lf", (*vector + i));
        }

        printf("Enter the matrix values row by row (space-separated):\n");
        for (int i = 0; i < (*row_count) * (*col_count); i++) {
            scanf("%lf", (*matrix + i));
        }
    } else {
        *vector = (double *)malloc((*col_count) * sizeof(double));
    }

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int rows_per_process = *row_count / comm_size;
    *local_matrix = (double *)malloc(rows_per_process * (*col_count) * sizeof(double));

    MPI_Bcast(row_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(col_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(*vector, *col_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(*matrix, rows_per_process * (*col_count), MPI_DOUBLE, *local_matrix, rows_per_process * (*col_count), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

/*
 * Function to perform matrix-vector multiplication.
 */
void BCS221050_Multiplication(int col_count, int rows_per_process, double *local_matrix, double *vector, double *local_result) {
    for (int i = 0; i < rows_per_process; i++) {
        local_result[i] = 0.0;
        for (int j = 0; j < col_count; j++) {
            local_result[i] += local_matrix[i * col_count + j] * vector[j];
        }
    }
}

/*
 * Function to gather and print the final result.
 */
void BCS221050_Print(int my_rank, double *local_result, int row_count, int comm_size) {
    double *global_result = NULL;
    int rows_per_process = row_count / comm_size;

    if (my_rank == 0) {
        global_result = (double *)malloc(row_count * sizeof(double));
    }

    MPI_Gather(local_result, rows_per_process, MPI_DOUBLE, global_result, rows_per_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("Resulting vector after matrix-vector multiplication:\n");
        for (int i = 0; i < row_count; i++) {
            printf("Result[%d]: %.2lf\n", i, global_result[i]);
        }
        free(global_result);
    }
}

/*
 * Function to calculate and display the maximum elapsed time across all processes.
 */
void BCS221050_Mytime(int my_rank, double start_time, double end_time) {
    double elapsed_time = end_time - start_time;
    double max_elapsed_time;

    MPI_Reduce(&elapsed_time, &max_elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("Maximum time taken by any process: %.6lf seconds\n", max_elapsed_time);
    }
}
