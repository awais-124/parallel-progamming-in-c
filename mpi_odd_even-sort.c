#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

const int RMAX = 100;

void Usage(char* program);
void Print_list(int local_A[], int local_n, int rank);
void Merge_low(int local_A[], int temp_B[], int temp_C[], int local_n);
void Merge_high(int local_A[], int temp_B[], int temp_C[], int local_n);
void Generate_list(int local_A[], int local_n, int my_rank);
int  Compare(const void* a_p, const void* b_p);
void Get_args(int argc, char* argv[], int* global_n_p, int* local_n_p, char* gi_p, int my_rank, int p, MPI_Comm comm);
void Sort(int local_A[], int local_n, int my_rank, int p, MPI_Comm comm);
void Odd_even_iter(int local_A[], int temp_B[], int temp_C[], int local_n, int phase, int even_partner, int odd_partner, int my_rank, int p, MPI_Comm comm);
void Print_local_lists(int local_A[], int local_n, int my_rank, int p, MPI_Comm comm);
void Print_global_list(int local_A[], int local_n, int my_rank, int p, MPI_Comm comm);
void Read_list(int local_A[], int local_n, int my_rank, int p, MPI_Comm comm);

int main(int argc, char* argv[]) {
   int my_rank, p;
   char g_i;
   int *local_A;
   int global_n;
   int local_n;
   MPI_Comm comm;

   MPI_Init(&argc, &argv);
   comm = MPI_COMM_WORLD;
   MPI_Comm_size(comm, &p);
   MPI_Comm_rank(comm, &my_rank);

   Get_args(argc, argv, &global_n, &local_n, &g_i, my_rank, p, comm);
   local_A = (int*) malloc(local_n*sizeof(int));
   if (g_i == 'g') {
      Generate_list(local_A, local_n, my_rank);
      Print_local_lists(local_A, local_n, my_rank, p, comm);
   } else {
      Read_list(local_A, local_n, my_rank, p, comm);
   }

   Sort(local_A, local_n, my_rank, p, comm);
   Print_global_list(local_A, local_n, my_rank, p, comm);
   free(local_A);
   MPI_Finalize();
   return 0;
}

void Generate_list(int local_A[], int local_n, int my_rank) {
   int i;
   srandom(my_rank+1);
   for (i = 0; i < local_n; i++)
      local_A[i] = random() % RMAX;
}

void Usage(char* program) {
   fprintf(stderr, "usage:  mpirun -np <p> %s <g|i> <global_n>\n", program);
   fprintf(stderr, "   - p: the number of processes \n");
   fprintf(stderr, "   - g: generate random, distributed list\n");
   fprintf(stderr, "   - i: user will input list on process 0\n");
   fprintf(stderr, "   - global_n: number of elements in global list");
   fprintf(stderr, " (must be evenly divisible by p)\n");
   fflush(stderr);
}

void Get_args(int argc, char* argv[], int* global_n_p, int* local_n_p, char* gi_p, int my_rank, int p, MPI_Comm comm) {
   if (my_rank == 0) {
      if (argc != 3) {
         Usage(argv[0]);
         *global_n_p = -1;
      } else {
         *gi_p = argv[1][0];
         if (*gi_p != 'g' && *gi_p != 'i') {
            Usage(argv[0]);
            *global_n_p = -1;
         } else {
            *global_n_p = atoi(argv[2]);
            if (*global_n_p % p != 0) {
               Usage(argv[0]);
               *global_n_p = -1;
            }
         }
      }
   }

   MPI_Bcast(gi_p, 1, MPI_CHAR, 0, comm);
   MPI_Bcast(global_n_p, 1, MPI_INT, 0, comm);

   if (*global_n_p <= 0) {
      MPI_Finalize();
      exit(-1);
   }
   *local_n_p = *global_n_p/p;
}

void Read_list(int local_A[], int local_n, int my_rank, int p, MPI_Comm comm) {
   int i;
   int *temp;
   if (my_rank == 0) {
      temp = (int*) malloc(p*local_n*sizeof(int));
      printf("Enter the elements of the list\n");
      for (i = 0; i < p*local_n; i++)
         scanf("%d", &temp[i]);
   } 
   MPI_Scatter(temp, local_n, MPI_INT, local_A, local_n, MPI_INT, 0, comm);
   if (my_rank == 0)
      free(temp);
}

void Print_global_list(int local_A[], int local_n, int my_rank, int p, MPI_Comm comm) {
   int* A;
   int i, n;
   if (my_rank == 0) {
      n = p*local_n;
      A = (int*) malloc(n*sizeof(int));
      MPI_Gather(local_A, local_n, MPI_INT, A, local_n, MPI_INT, 0, comm);
      printf("Global list:\n");
      for (i = 0; i < n; i++)
         printf("%d ", A[i]);
      printf("\n\n");
      free(A);
   } else {
      MPI_Gather(local_A, local_n, MPI_INT, A, local_n, MPI_INT, 0, comm);
   }
}

int Compare(const void* a_p, const void* b_p) {
   int a = *((int*)a_p);
   int b = *((int*)b_p);
   if (a < b)
      return -1;
   else if (a == b)
      return 0;
   else
      return 1;
}
