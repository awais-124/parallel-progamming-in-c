/* File:     mpi_trap1.c
 * Purpose:  Use MPI to implement a parallel version of the trapezoidal 
 *           rule.  In this version the endpoints of the interval and
 *           the number of trapezoids are hardwired.
 *
 * Input:    None.
 * Output:   Estimate of the integral from a to b of f(x)
 *           using the trapezoidal rule and n trapezoids.
 *
 * Compile:  mpicc -g -Wall -o mpi_trap1 mpi_trap1.c
 * Run:      mpiexec -n <number of processes> ./mpi_trap1
 *
 * Algorithm:
 *    1.  Each process calculates "its" interval of
 *        integration.
 *    2.  Each process estimates the integral of f(x)
 *        over its interval using the trapezoidal rule.
 *    3a. Each process != 0 sends its integral to 0.
 *    3b. Process 0 sums the calculations received from
 *        the individual processes and prints the result.
 *
 * Note:  f(x), a, b, and n are all hardwired.
 *
 * IPP:   Section 3.2.2 (pp. 96 and ff.)
 */
#include <stdio.h>

/* We'll be using MPI routines, definitions, etc. */
#include <mpi.h>

/* Calculate local integral  */
double Trap(double left_endpt, double right_endpt, int trap_count, 
   double base_len);    

/* Function we're integrating */
double f(double x); 

void BCS221057_Print(double total_int, int n, double a, double b);

void BCS221057_Input(int my_rank, int comm_sz, int* n, double *a, double *b);



int main(void) {
   int my_rank, comm_sz, n, local_n;   
   double a, b, h, local_a, local_b;
   double local_int, total_int;
   int source; 

   /* Let the system do what it needs to start up MPI */
   MPI_Init(NULL, NULL);

   /* Get my process rank */
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   /* Find out how many processes are being used */
   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
   
   
   
   BCS221057_Input(my_rank, comm_sz, &n, &a, &b);
   
   

   h = (b-a)/n;          /* h is the same for all processes */
   local_n = n/comm_sz;  /* So is the number of trapezoids  */

   /* Length of each process' interval of
    * integration = local_n*h.  So my interval
    * starts at: */
   local_a = a + my_rank*local_n*h;
   local_b = local_a + local_n*h;
   local_int = Trap(local_a, local_b, local_n, h);

   /* Add up the integrals calculated by each process */
   if (my_rank != 0) { 
      MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, 
            MPI_COMM_WORLD);
   } else {
      total_int = local_int;
      printf("\nThe trapezoid area calculated by Process %d of %d is %lf", my_rank,comm_sz,local_int);
      for (source = 1; source < comm_sz; source++) {
         MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         total_int += local_int;
         printf("\nThe trapezoid area calculated by Process %d of %d is %lf", my_rank,comm_sz,local_int);
      }
   } 

   /* Print the result */
   if (my_rank == 0) {
   
      BCS221057_Print(total_int, n, a, b);
      
   }

   /* Shut down MPI */
   MPI_Finalize();

   return 0;
} /*  main  */

/*------------------------------------------------------------------
 * Function:     Trap
 * Purpose:      Serial function for estimating a definite integral 
 *               using the trapezoidal rule
 * Input args:   left_endpt
 *               right_endpt
 *               trap_count 
 *               base_len
 * Return val:   Trapezoidal rule estimate of integral from
 *               left_endpt to right_endpt using trap_count
 *               trapezoids
 */
double Trap(
      double left_endpt  /* in */, 
      double right_endpt /* in */, 
      int    trap_count  /* in */, 
      double base_len    /* in */) {
   double estimate, x; 
   int i;

   estimate = (f(left_endpt) + f(right_endpt))/2.0;
   for (i = 1; i <= trap_count-1; i++) {
      x = left_endpt + i*base_len;
      estimate += f(x);
   }
   estimate = estimate*base_len;

   return estimate;
} /*  Trap  */


/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 */
double f(double x) {
   return x*x;
} /* f */





void BCS221057_Print(double total_int, int n, double a, double b) {

      printf("\nWith n = %d trapezoids, our estimate\n", n);
      printf("of the integral from %f to %f = %.15e\n",
          a, b, total_int);
}

void BCS221057_Input(int my_rank, int comm_sz, int* n, double *a, double *b) {

    if (my_rank == 0) {
    
        printf("Enter values of a, b, and n (one by one):  ");
        scanf("%lf %lf %d", a, b, n);
        for (int dest = 1; dest < comm_sz; dest++) {
            MPI_Send(a, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            MPI_Send(b, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            MPI_Send(n, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
            
        }
    } else 
    {
        MPI_Recv(a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}
