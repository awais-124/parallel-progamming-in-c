#include <stdio.h>
#include <mpi.h>

double Trap(double left_endpt, double right_endpt, int trap_count, double base_len);
double f(double x);

void BCS221098_Input(double* a, double* b, int* n, int my_rank, int comm_sz);
void BCS221098_Print(double total_int);


int main(void) {
    int my_rank, comm_sz, n, local_n;
    double local_int, total_int;
    double a, b, h, local_a, local_b;

    MPI_Init(NULL, NULL); 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);	
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    /* PART - A */
    BCS221098_Input(&a, &b, &n, my_rank, comm_sz);	

    h = (b - a) / n;
    local_n = n / comm_sz; 
    local_a = a + my_rank * local_n * h;	
    local_b = local_a + local_n * h;	
    local_int = Trap(local_a, local_b, local_n, h);	

    if (my_rank != 0) {
        MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else { 
        printf("Process %d of %d: My Area  is %lf\n",my_rank, comm_sz, local_int);
        for (int process_rank = 1; process_rank < comm_sz; process_rank++) {
            MPI_Recv(&local_int, 1, MPI_DOUBLE, process_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_int += local_int;
            printf("Process %d of %d: My Area is is %lf\n", process_rank,  comm_sz, local_int);
        }
        /* sum all local integrals */
        total_int+=local_int;
    }

    /* PART - C */
    if (my_rank == 0) { BCS221098_Print(total_int);}
    MPI_Finalize();
    return 0;
}


void BCS221098_Input(double* a, double* b, int* n, int my_rank, int comm_sz) {
    if (my_rank != 0) { 	
        MPI_Recv(n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    } 

    if( my_rank == 0 ) {
        printf("Enter a  b and n in sequence: ");
        scanf("%lf %lf %d", a, b, n);
        
        for (int process_rank= 1, process_rank < comm_sz, process_rank++) {
            MPI_Send(n, 1, MPI_INT, process_rank, 0, MPI_COMM_WORLD);
            MPI_Send(b, 1, MPI_DOUBLE, process_rank, 0, MPI_COMM_WORLD);
            MPI_Send(a, 1, MPI_DOUBLE, process_rank, 0, MPI_COMM_WORLD);
        }
    }
}

void BCS221098_Print(double total_int) {
    printf("\nTOTAL AREA IS  =  %lf", total_int);
}

double Trap(double left_endpt, double right_endpt, int trap_count, double base_len) {
    double integral, x;
    int i;
    integral = (f(left_endpt) + f(right_endpt)) / 2.0;
    for (i = 1; i <= trap_count - 1; i++) {
        x = left_endpt + i * base_len;
        integral += f(x);
    }
    integral = integral *  base_len;
    return integral;
}

double f(double x) {
    return x * x;
}
