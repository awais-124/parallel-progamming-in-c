/* BCS221093 - Muhammad Awais
   PDC Assignment 03 - (26/12/24) */
#include <stdio.h>
#include <mpi.h>

/*FUNCTION PROTOTYPES*/
void BCS221093_Input(double* a, double* b, int* n, int my_rank, int comm_sz);
void BCS221093_Print(double total_int);
double Trap(double left_endpt, double right_endpt, int trap_count, double base_len);
double f(double x);


/*MAIN FUNCTION*/
int main(void) {
    int my_rank, comm_sz, n, local_n;
    double a, b, h, local_a, local_b;
    double local_int, total_int;

    MPI_Init(NULL, NULL); /* allocate necessary resources to this program */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);	/* initialize rank of current process */
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);	/* initialize number of processes, passes via terminal */

	
	/* Part A : P0 takes inputs and send them to other Processes using MPI_Send,
	 Other Processes receive inputs using MPI_Recv */
    BCS221093_Input(&a, &b, &n, my_rank, comm_sz);	

    h = (b - a) / n; /* width of one trapeziod */
    local_n = n / comm_sz; /* number of trapeziods for a specific process, equal for all*/


    local_a = a + my_rank * local_n * h;	/* starting point for specific process. */
    local_b = local_a + local_n * h;		/* ending point for specific process. */

    
    /* function to calculate the area(integral) for specific process(my_rank)*/
    local_int = Trap(local_a, local_b, local_n, h);	

	
    /* creating sender-receiver split. Process 0 will receive local_integrals 
       from all other processes and add them to calculate final area */
    if (my_rank != 0) {
    	/* send your integral to P0(master core) */
        MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else { /* rank == 0 */
        total_int = local_int;
        printf("PROCESS 0 OF %d : TRAPEZIOD AREA CALCULATED IS (%.10f)\n", comm_sz, local_int);

        for (int source_rank = 1; source_rank < comm_sz; source_rank++) {
            MPI_Recv(&local_int, 1, MPI_DOUBLE, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_int += local_int;
            printf("PROCESS %d OF %d : TRAPEZIOD AREA CALCULATED IS (%.10f)\n", source_rank, comm_sz, local_int);
        }
		
		/* Print total result(area),
		   Only Process 0 will print, 
		   this line is inside else block */
        BCS221093_Print(total_int);
    }

    /* Release resources */
    MPI_Finalize();

    return 0;
}

/* FUNCTION DEFINITIONS */
void BCS221093_Input(double* a, double* b, int* n, int my_rank, int comm_sz) {
    if (my_rank == 0) { 	/* Master core functionality */
    
        printf("Enter the value for a(starting point): ");
        scanf("%lf", a);
        
        printf("Enter the value for b(ending point): ");
        scanf("%lf", b);
        
        printf("Enter the value for n(number of trapeziods)(n > 0): ");
        /* n must be positive as it is a count, adding check for that */
        do {
        	scanf("%d", n);
        	if (n <= 0) printf("Plz enter n > 0: ");
        } while(n <= 0); // if n is 0 or negative, take input again
        
        
        /* Need to send inputs to all other processes, using loop for that
           Sending each input separately using MPI_Send() */
        for (int dest_rank = 1; dest_rank < comm_sz; dest_rank++) {
            MPI_Send(a, 1, MPI_DOUBLE, dest_rank, 0, MPI_COMM_WORLD);
            MPI_Send(b, 1, MPI_DOUBLE, dest_rank, 0, MPI_COMM_WORLD);
            MPI_Send(n, 1, MPI_INT, dest_rank, 0, MPI_COMM_WORLD);
        }
        
    } else {	/* Worker cores receiving inputs from P0(master)*/
        MPI_Recv(a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

/* function printin final result */
void BCS221093_Print(double total_int) {
    printf("PROCESS 0 SAYS \"TOTAL AREA CALCULATED =  %.10f\"\n", total_int);
}

/* funtion to calculate area, copied from mpi_trap1.c file*/
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

/* square function, copied from mpi_trap1.c file */
double f(double x) {
    return x * x;
}
