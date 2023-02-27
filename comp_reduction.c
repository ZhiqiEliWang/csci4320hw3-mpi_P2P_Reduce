#include <stdio.h>
#include <mpi.h>



MPI_P2P_Reduce(void* send_data, // each process's partition of task array
    void* recv_data, // the result of reduction, for root only
    int count, // length of the send_data
    MPI_Datatype datatype,
    int root,
    MPI_Comm communicator){ // this function takes the same inputs as MPI_reduce, except MPI_Op is set to MPI_SUM


    // -----------1. Each rank computes sum over local data array.---------------
    for (int i=0; i<)


}


int main(int argc, char** argv){
  // Initialize the MPI environment
    int world_rank, world_size; // init world rank and size
    MPI_Init(&argc, &argv);

    // Find out rank, size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // size of a array is determined by how many nodes are working on this task
    arrSize = 2**30 / world_rank; 

    MPI_LONG_LONG* bigArr = malloc(sizeof(MPI_LONG_LONG)*world_size);

    for (int i=0; i<arrSize; i++){
        bigArr[i] = world_rank * arrSize;
    }


    // calling MPI_P2P_Reduce
    MPI_LONG_LONG global_sum; // only 0th rank process will have this value modified 
    MPI_P2P_Reduce(bigArr, &global_sum, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    free(bigArr);
    MPI_Finalize();
    return 0;
}