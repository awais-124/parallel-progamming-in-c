#include <stdio.h>
#include <mpi.h>

double Trap(double left_endpt, double right_endpt, int trap_count, double base_len);
double f(double x);
void BCS221099_Input(double* a, double* b, int* n, int my_rank, int comm_sz);
void BCS221099_Print(double total_int);


int main(void) {
    int my_rank, comm_sz, n, local_n;
    double local_int, total_int;
    double a, b, h, local_a, local_b;

    MPI_Init(NULL, NULL); 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);	
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    /* PART A Code */
    BCS221099_Input(&a, &b, &n, my_rank, comm_sz);	

    h = (b - a) / n;
    local_n = n / comm_sz; 
    local_a = a + my_rank * local_n * h;	
    local_b = local_a + local_n * h;	
    local_int = Trap(local_a, local_b, local_n, h);	

    if (my_rank != 0) {
        MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else { 
        printf("Local integral of Process %d of %d is %lf\n",my_rank, comm_sz, local_int );
        for (int q = 1; q < comm_sz; q++) {
            MPI_Recv(&local_int, 1, MPI_DOUBLE, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_int += local_int;
            printf("Local integral of Process %d of %d is %lf\n", q, comm_sz, local_int);
        }
        total_int+=local_int;
    }


    // PART C Code
    if (my_rank == 0){
        BCS221099_Print(total_int);
    }

    MPI_Finalize();

    print("\n");
    return 0;
}


void BCS221093_Input(double* a, double* b, int* n, int my_rank, int comm_sz) {
    if (my_rank == 0) { 	
        
        printf("Enter a and b: ");
        scanf("%lf %lf", a, b);

        printf("Enter value for n: ");
        scanf("%d", n);
        
        for (int i= 1, i < comm_sz, i++) {
            MPI_Send(n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(b, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send(a, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        
    } else {
        MPI_Recv(n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void BCS221099_Print(double total_int) {
    printf("\nTotal integral is =  %lf", total_int);
}

double f(double x) {
    return x * x;
}

double Trap(double left_endpt, double right_endpt, int trap_count, double base_len) {
    double approx, x;
    int i;

    approx = (f(left_endpt) + f(right_endpt)) / 2.0;
    for (i = 1; i <= trap_count - 1; i++) {
        x = left_endpt + i * base_len;
        approx += f(x);
    }
    approx *=  base_len;
    return approx;
}

