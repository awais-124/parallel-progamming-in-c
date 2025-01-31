#include <stdio.h>
#include <mpi.h>

void BCS221093_Input(double* a, double* b, int* n, int my_rank, int comm_sz);
void BCS221093_Print(double total_int);
double Trap(double left_endpt, double right_endpt, int trap_count, double base_len);
double f(double x);

int main(void) {
    int my_rank, comm_sz, n, local_n;
    double a, b, h, local_a, local_b;
    double local_int, total_int;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    BCS221093_Input(&a, &b, &n, my_rank, comm_sz);

    h = (b - a) / n;
    local_n = n / comm_sz;

    local_a = a + my_rank * local_n * h;
    local_b = local_a + local_n * h;

    local_int = Trap(local_a, local_b, local_n, h);

    if (my_rank != 0) {
        MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else {
        total_int = local_int;
        printf("The trapezoid area calculated by Process 0 of %d is %.15e\n", comm_sz, local_int);

        for (int source = 1; source < comm_sz; source++) {
            MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_int += local_int;
            printf("The trapezoid area received by Process %d of %d is %.15e\n", source, comm_sz, local_int);
        }

        BCS221093_Print(total_int);
    }

    MPI_Finalize();

    return 0;
}

void BCS221093_Input(double* a, double* b, int* n, int my_rank, int comm_sz) {
    if (my_rank == 0) {
        printf("Enter the values for a, b, and n (e.g., 0 36 12): ");
        scanf("%lf %lf %d", a, b, n);
        for (int dest = 1; dest < comm_sz; dest++) {
            MPI_Send(a, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            MPI_Send(b, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            MPI_Send(n, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void BCS221093_Print(double total_int) {
    printf("The total integral is %.15e\n", total_int);
}

double Trap(double left_endpt, double right_endpt, int trap_count, double base_len) {
    double estimate, x;
    int i;

    estimate = (f(left_endpt) + f(right_endpt)) / 2.0;
    for (i = 1; i <= trap_count - 1; i++) {
        x = left_endpt + i * base_len;
        estimate += f(x);
    }
    estimate = estimate * base_len;

    return estimate;
}

double f(double x) {
    return x * x;
}
