#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>



// this function takes the same inputs as MPI_reduce, except MPI_Op is set to MPI_SUM
int MPI_P2P_Reduce(long long int* send_data, // each process's partition of task array
    long long int* recv_data, // the result of reduction, for root only
    int count, // length of the send_data
    MPI_Datatype datatype, // MPI_LONG_LONG in our case
    int root, // 0 in our case
    MPI_Comm communicator) // we only have one communicator
    {

    long long int  local_sum = 0;
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
        if ((self_rank / stride) % 2 == 0 && self_rank % stride == 0){ // even ranks after stride: receiver
            long long int recv_buf;
            MPI_Request recv_req;
            MPI_Status recv_status;
            MPI_Irecv(&recv_buf, 1, MPI_LONG_LONG, self_rank+1, MPI_ANY_TAG, communicator, &recv_req);
            MPI_Wait(&recv_req , &recv_status);
            local_sum += recv_buf; // perform pairwise sum here
        }
        else if (self_rank % stride == 0){ // odd ranks after stride: sender
            MPI_Request send_req;
            MPI_Status send_status;
            MPI_Isend(&local_sum, 1, MPI_LONG_LONG, self_rank-1, 0, communicator, &send_req);
            MPI_Wait(&send_req , &send_status);
        }
        stride *= 2;
        printf("rank: %d hits barrier\n", self_rank);
        MPI_Barrier(communicator); // sync here
    }
    return 0;
}


int main(int argc, char** argv){
  // Initialize the MPI environment
    int world_rank, world_size; // init world rank and size
    MPI_Init(&argc, &argv);

    // Find out rank, size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // size of a array is determined by how many nodes are working on this task
    int arrSize = 2<<30 / world_rank; 

    long long int* bigArr = malloc(sizeof(long long int)*world_size);

    for (int i=0; i<arrSize; i++){
        bigArr[i] = world_rank * arrSize;
    }

    // calling MPI_P2P_Reduce
    long long int global_sum; // only 0th rank process will have this value modified 
    MPI_P2P_Reduce(bigArr, &global_sum, arrSize, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ // only root has the result and it prints it here
        long long int correct_answer = 576460751766552576;
        if (correct_answer - global_sum == 0){
            printf("result is correct: %lld\n", correct_answer);
        }
        else{
            printf("correct answer: %lld\nmy answer:      %lld\n",correct_answer, global_sum);
        }
    }

    free(bigArr);
    MPI_Finalize();
    return 0;
}