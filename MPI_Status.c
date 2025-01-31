#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    int rank, size, data[10];
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {  // Sender
        int send_data[5] = {1, 2, 3, 4, 5};
        MPI_Send(send_data, 5, MPI_INT, 1, 99, MPI_COMM_WORLD);
    } 
    else if (rank == 1) {  // Receiver (unknown sender and size)
        MPI_Recv(data, 10, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        // Get the actual count of received elements
        int count;
        MPI_Get_count(&status, MPI_INT, &count);

        printf("Received %d elements from process %d with tag %d\n", 
               count, status.MPI_SOURCE, status.MPI_TAG);
    }

    MPI_Finalize();
    return 0;
}
