/* BCS221093 - Muhammad Awais
   PDC Assignment 04 - (15 / 01 / 25) */

/*
    USED mpi_mat_vect_mult.c and mpi_mat_vect_time.c 
    files from ippsource codes for reference.
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/* PROTOTYPES OF REQUIRED FUNCTIONS */
void BCS221093_Input(int, int *, int *, double **, double **, double **, int);
void BCS221093_Multiplication(int, int, int, int, double *, double *, double *);
void BCS221093_Print(int, int, int, double *);
void BCS221093_Mytime(int, double, double);

/* MAIN FUNCTION */
int main() {
    int my_rank, comm_size;
    int rows=0, cols=0;

    /* a is matrix and x is vector*/
    double *a = NULL, *x = NULL;

    /* local_a refers to part of matrix local to process and same for local_y*/
    double *local_a = NULL, *local_y = NULL;

    double start_time, end_time;

    MPI_Init(NULL, NULL);   /* allocate necessary resources to this program */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /* initialize rank of current process */
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);  /* initialize number of processes, passes via terminal */

    /* THIS CHECK IS ADDED TO MAKE SURE AT LEAST 2 PROCESSES ARE BEING USED
    SO THAT PARALLELISM CAN BE ACHIEVED. ALSO USING MPI_BCAST AND MPI_GATHER WITH
    ONE PROCESS CAN CAUSE ERRORS. TO HANDLE SUCH ERRORS, BELOW CHECK IS ADDED.
    MPI_ABORT() ABORTS THE PROGRAM WITH PROPER TERMINATION */
    if (my_rank == 0 && comm_size < 2) {
        printf("RUN THE PROGRAM WITH AT LEAST TWO PROCESSES.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    /* ADDING BARRIERS SO THAT ALL PROGRAM STARTS CALCULATING
    TIME AT THE SAME TIME */
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime(); // start time calculation for all processes

    /* CALL INPUT FUNCTION */
    BCS221093_Input(my_rank, &rows, &cols, &a, &x, &local_a, comm_size); // taking input and scattering data

    int rows_per_process = rows / comm_size;
    local_y = (double *)malloc(rows_per_process * sizeof(double));

    /* MATRIX-VECTOR MULTIPLICATION FUNCTION */
    BCS221093_Multiplication(my_rank, comm_size, rows_per_process, cols, local_a, x, local_y);

    /* FUNCTION TO PRINT RESULTS */
    BCS221093_Print(my_rank, comm_size, rows, local_y);

    end_time = MPI_Wtime();

    /* FUNCTION TO CALCULATE AND PRINT EXECUTION TIME */
    BCS221093_Mytime(my_rank, start_time, end_time);

    /* RELEASE THE DYNAMICALLY ALLOCATED MEMORY */
    if (my_rank == 0) {
        /* only process 0 assigned memory to a, thats
        why only P0 cleans it */
        free(a);
    }
    free(x);
    free(local_a);
    free(local_y);

    MPI_Finalize();
    return 0;
}
/* MAIN ENDED */

/* FUNCTION DEFINITIONS */

void BCS221093_Input(int my_rank, int *rows, int *cols, double **a, double **x, double **local_a, int comm_size) {

    /* MEMORY ALLOCATION FOR VECTOR */
    *x = (double *)malloc((*cols) * sizeof(double));

    if (my_rank == 0) {
        /* INPUT ROWS AND COLUMNS (ORDER OF MATRIX)
            KEEP TAKING INPUT UNTILL POSITIVE NO IS ENTERED */
        while (*rows <= 0) {
            printf("ENTER NUMBER OF ROWS ( > 0 ): ");
            scanf("%d", rows);
            if (*rows <= 0) {
                printf("INVALID INPUT. ENTER A POSITIVE NUMBER PLEASE...\n");
            }
        }

        while (*cols <= 0) {
            printf("ENTER NUMBER OF COLUMNS ( > 0 ): ");
            scanf("%d", cols);
            if (*cols <= 0) {
                printf("INVALID INPUT. ENTER A POSITIVE NUMBER PLEASE...\n");
            }
        }


        /* MEMORY ALLOCAION FOR MATRIX */
        *a = (double *)malloc((*rows) * (*cols) * sizeof(double));

        /* TAKING INPUT MATRIX VALUES */
        printf("ENTER THE MATRIX VALUES (%dx%d):\n", *rows, *cols);
        for (int i = 0; i < *rows; i++) {
            for (int j = 0; j < *cols; j++) {
                printf("VALUE AT [%d][%d]: ", i, j);
                scanf("%lf", (*a + i * (*cols) + j));
            }
        }

        /* TAKING INPUT VECTOR VALUES */
        printf("ENTER THE VECTOR VALUES (%dx1):\n", *cols);
        for (int i = 0; i < *cols; i++) {
            printf("VALUE AT [%d]: ", i);
            scanf("%lf", (*x + i));
        }
    }

    /* USING MPI_Bcast() to share ROWS, COLUMNS, AND VECTOR with ALL PROCESSES */
    MPI_Bcast(rows, 1, MPI_INT, 0, MPI_COMM_WORLD); /* ROWS */
    MPI_Bcast(cols, 1, MPI_INT, 0, MPI_COMM_WORLD); /* COLUMNS */
    MPI_Bcast(*x, *cols, MPI_DOUBLE, 0, MPI_COMM_WORLD); /* VECTORS */

    /* ALLOCATE MEMORY FOR LOCAL PART OF MATRIX */
    int rows_per_process = (*rows) / comm_size;
    *local_a = (double *)malloc(rows_per_process * (*cols) * sizeof(double));

    /* SCATTERING MATRIX VALUES EACH PROCESS GETS (rows_per_process X cols) values */
    MPI_Scatter(*a, rows_per_process * (*cols), MPI_DOUBLE, *local_a, rows_per_process * (*cols), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void BCS221093_Multiplication(int my_rank, int comm_size, int rows_per_process, int cols, double *local_a, double *x, double *local_y) {
    for (int i = 0; i < rows_per_process; i++) {
        local_y[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            local_y[i] += local_a[i * cols + j] * x[j];
        }
    }
}

void BCS221093_Print(int my_rank, int comm_size, int rows, double *local_y) {
    double *global_y = NULL;
    int rows_per_process = rows / comm_size;

    if (my_rank == 0) {
        global_y = (double *)malloc(rows * sizeof(double));
    }

    /* COLLECT LOCAL RESULTS INTO GLOBAL VECTOR Y */
    MPI_Gather(local_y, rows_per_process, MPI_DOUBLE, global_y, rows_per_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        /* PRINTING RESULTANT RESULT VECTOR */
        printf("RESULTANT VECTOR y:\n");
        for (int i = 0; i < rows; i++) {
            printf("y[%d]: %.2lf\n", i, global_y[i]);
        }
        free(global_y);
    }
}

void BCS221093_Mytime(int my_rank, double start_time, double end_time) {
    double local_time = end_time - start_time;
    double max_time;

    /* MPI_REDUCE TO FIND THE MAXIMUM TIME AMONG ALL PROCESSES */
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("MAXIMUM TIME TAKEN BY ANY PROCESS: %.6lf SECONDS\n", max_time);
    }
}
