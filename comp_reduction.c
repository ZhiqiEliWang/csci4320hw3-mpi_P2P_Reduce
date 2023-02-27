#include <stdio.h>
#include <mpi.h>



// this function takes the same inputs as MPI_reduce, except MPI_Op is set to MPI_SUM
MPI_P2P_Reduce(void* send_data, // each process's partition of task array
    void* recv_data, // the result of reduction, for root only
    int count, // length of the send_data
    MPI_Datatype datatype, // MPI_LONG_LONG in our case
    int root, // 0 in our case
    MPI_Comm communicator) // we only have one communicator
    {

    MPI_LONG_LONG local_sum = 0;
    int self_rank;
    MPI_Comm_rank(communicator, &self_rank); // get current rank number
    int comm_size;
    MPI_Comm_size(communicator, &comm_size); // get total number of nodes

    // -----------1. Each rank computes sum over local data array.---------------
    for (int i=0; i<count; i++){
        local_sum += send_data[i];
    }

    // --------------2. Compute pairwise sums between MPI ranks-------------------
    int stride = 1;
    while (stride < comm_size){
        if ((self_rank / stride) % 2){ // even ranks after stride: receiver
            MPI_LONG_LONG recv_buf;
            MPI_Request recv_req;
            MPI_Status recv_status;
            MPI_Irecv(&recv_buf, 1, MPI_LONG_LONG, self_rank+1, MPI_ANY_TAG, communicator, &recv_req);
            MPI_Wait(&recv_buf , &recv_status);
            local_sum += recv_buf; // perform pairwise sum here
        }
        else{ // odd ranks after stride: sender
            MPI_Request send_req;
            MPI_Status send_status;
            MPI_Isend(&local_sum, 1, MPI_LONG_LONG, self_rank-1, MPI_ANY_TAG, communicator, &send_req);
            MPI_Wait(&send_buf , &send_status);
        }
        stride *= 2;
        MPI_Barrier(communicator); // sync here
    }
    return;
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
    MPI_P2P_Reduce(bigArr, &global_sum, arrSize, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ // only root has the result and it prints it here
        MPI_LONG_LONG currect_answer = 576460751766552576;
        if (currect_answer - global_sum == 0){
            printf("result is correct: %lld\n", currect_answer);
        }
        else{
            printf("correct answer: %lld\nmy answer:      %lld\n");
        }
    }

    free(bigArr);
    MPI_Finalize();
    return 0;
}