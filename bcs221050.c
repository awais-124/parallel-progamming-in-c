/* File:     mpi_trap1.c
 * Purpose:  Use MPI to implement a parallel version of the trapezoidal 
 *           rule. In this version, the endpoints of the interval and
 *           the number of trapezoids are predefined.
 *
 * Input:    None.
 * Output:   Estimate of the integral from start to end of f(x)
 *           using the trapezoidal rule with n trapezoids.
 *
 * Compile:  mpicc -g -Wall -o mpi_trap1 mpi_trap1.c
 * Run:      mpiexec -n <number of processes> ./mpi_trap1
 *
 * Algorithm:
 *    1.  Each process determines its interval of
 *        integration.
 *    2.  Each process computes the integral of f(x)
 *        over its interval using the trapezoidal rule.
 *    3a. Each process (except process 0) sends its integral to process 0.
 *    3b. Process 0 aggregates the results received from
 *        all processes and displays the final result.
 *
 * Note:  f(x), start, end, and n are predefined.
 */
#include <stdio.h>

/* Include MPI library */
#include <mpi.h>

// Function Prototypes
void BCS221050_Input(double *start, double *end, int *num_trapezoids, int process_rank, int total_processes);
void BCS221050_Print(double total_integral, int num_trapezoids, double start, double end);

/* Function to calculate local integral */
double ComputeLocalIntegral(double interval_start, double interval_end, int trapezoid_count, 
   double interval_length);

/* Function being integrated */
double Function(double x);

// Function Definitions
#include <stdio.h>

void BCS221050_Input(double *start, double *end, int *num_trapezoids, int process_rank, int total_processes) {
 if (process_rank == 0) {
	 // Prompt user to input the interval and number of trapezoids
	 printf("Enter the start of the interval: ");
	 scanf("%lf", start);
	 
	 printf("Enter the end of the interval: ");
	 scanf("%lf", end);
	 
	 printf("Enter the number of trapezoids: ");
	 scanf("%d", num_trapezoids);
      // Broadcast input values to all other processes
      for (int dest = 1; dest < total_processes; dest++) {
         MPI_Send(start, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);

         MPI_Send(end, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
         MPI_Send(num_trapezoids, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
      }
   } else {
      // Other processes receive the input values from process 0
      MPI_Recv(start, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(end, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(num_trapezoids, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }  
   
}

void BCS221050_Print(double total_integral, int num_trapezoids, double start, double end) {
   // Display the computed integral result
   printf("\nUsing n = %d trapezoids, the estimated integral\n", num_trapezoids);
   printf("from %.2f to %.2f is %.15e\n", start, end, total_integral);
}

int main(void) {
   int process_rank, total_processes = 4, num_trapezoids, local_trapezoids;   
   double start, end, interval_length, local_start, local_end;
   double local_integral, total_integral;
   int source_process; 

   /* Initialize MPI environment */
   MPI_Init(NULL, NULL);

   /* Get rank of the current process */
   MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

   /* Get total number of processes */
   MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
   
   // Process 0 takes input from the user and share to other processes
   BCS221050_Input(&start, &end, &num_trapezoids, process_rank, total_processes);

   interval_length = (end - start) / num_trapezoids;          /* Step size */
   local_trapezoids = num_trapezoids / total_processes;  /* Number of trapezoids per process */

   /* Determine the interval for each process */
   local_start = start + process_rank * local_trapezoids * interval_length;
   local_end = local_start + local_trapezoids * interval_length;
   local_integral = ComputeLocalIntegral(local_start, local_end, local_trapezoids, interval_length);

   /* Gather and sum the integrals from all processes */
   if (process_rank != 0) { 
      MPI_Send(&local_integral, 1, MPI_DOUBLE, 0, 0, 
            MPI_COMM_WORLD); 
   } else {
      total_integral = local_integral;
      // Process 0 displays intermediate results as they are received
      printf("The local integral computed by Process %d of %d is %.15e\n", process_rank, total_processes, local_integral);

      for (source_process = 1; source_process < total_processes; source_process++) {
         MPI_Recv(&local_integral, 1, MPI_DOUBLE, source_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         total_integral += local_integral;

         // Display intermediate results from other processes
         printf("The local integral computed by Process %d of %d is %.15e\n", source_process, total_processes, local_integral);
      }
   } 

/* Display the final result */
   if (process_rank == 0) {
      BCS221050_Print(total_integral, num_trapezoids, start, end);
   }
   

   /* Finalize MPI environment */
   MPI_Finalize();

   return 0;
} /*  main  */


/*------------------------------------------------------------------
 * Function:     ComputeLocalIntegral
 * Purpose:      Calculate the integral for a subinterval using the trapezoidal rule.
 * Input args:   interval_start
 *               interval_end
 *               trapezoid_count 
 *               interval_length
 * Return val:   Trapezoidal rule estimate for the subinterval
 */
double ComputeLocalIntegral(
      double interval_start  /* in */, 
      double interval_end /* in */, 
      int    trapezoid_count  /* in */, 
      double interval_length    /* in */) {
   double estimate, x; 
   int i;

   estimate = (Function(interval_start) + Function(interval_end)) / 2.0;
   for (i = 1; i <= trapezoid_count - 1; i++) {
      x = interval_start + i * interval_length;
      estimate += Function(x);
   }
   estimate = estimate * interval_length;

   return estimate;
} /*  ComputeLocalIntegral  */


/*------------------------------------------------------------------
 * Function:    Function
 * Purpose:     Define the function to be integrated
 * Input args:  x
 */
double Function(double x) {
   return x * x;
} /* Function */
