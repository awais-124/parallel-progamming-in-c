#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Send from one process to another (p2p Communication) @PARAMS: buffer for data, 
    count(number of elements), datatype, dest-rank,MPI tag for extra info, communicator */
    MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);

    /* Receive data from another process (p2p Communication) @PARAMS: buffer for data, 
    count (number of elements), datatype, source rank, MPI tag, communicator, status */
    MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, 
             MPI_Comm comm, MPI_Status *status);

    /* Broadcast data from one process to all processes (Collective Communication)
       @PARAMS: buffer (holds data to be broadcasted), count (number of elements),
       datatype, root process (source of broadcast), communicator */
    MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

    /* Scatter data from one process to all processes (Distributes data)
       @PARAMS: send buffer (data at root), send count (per process), send type,
       receive buffer (each process gets a portion), receive count, receive type,
       root process (who is scattering), communicator */
    MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                MPI_Datatype recvtype, int root, MPI_Comm comm);

    /* Gather data from all processes to one root process (Collective Communication)
       @PARAMS: send buffer (each process sends its data), send count, send type,
       receive buffer (only on root), receive count, receive type,
       root process (who gathers), communicator */
    MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, 
               int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);

    /* Performs Operations (e.g., sum, max, min)and stores the result in the root process
       @PARAMS: send buffer (data from each process), receive buffer (result at root),
       count (number of elements), datatype, operation (MPI_SUM, MPI_MAX, etc.),
       root process (who collects result), communicator */
    MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, 
               MPI_Op op, int root, MPI_Comm comm);

    /* Same as reduce, but all Ps gets result */
    MPI_Allreduce(&sum, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}
