#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "clockcycle.h"


// this function takes the same inputs as MPI_reduce, except MPI_Op is set to MPI_SUM
int MPI_P2P_Reduce(long long int* send_data, // each process's partition of task array
    long long int* recv_data, // the result of reduction, for root only
    int count, // length of the send_data
    MPI_Datatype datatype, // MPI_LONG_LONG in our case
    int root, // 0 in our case
    MPI_Comm communicator) // we only have one communicator
    {

    long long int local_sum = 0;
    int self_rank;
    MPI_Comm_rank(communicator, &self_rank); // get current rank number
    int comm_size;
    MPI_Comm_size(communicator, &comm_size); // get total number of nodes

    // -----------1. Each rank computes sum over local data array.---------------
    
    for (int i=0; i<count; i++){
        local_sum += send_data[i];
    }
    // printf("rank: %d, right before the barrier\n", self_rank);
    MPI_Barrier(communicator); // sync here

    // --------------2. Compute pairwise sums between MPI ranks-------------------
    int stride = 1;
    
    MPI_Request recv_req;
    MPI_Status recv_status;
    MPI_Request send_req;
    MPI_Status send_status;

    while (stride < comm_size){
        if (self_rank == 0){printf("---------stride: %d-----------\n", stride);}
        if ((self_rank / stride) % 2 == 0){ // even ranks after stride: receiver
            long long int recv_buf;
            MPI_Irecv(&recv_buf, 1, MPI_LONG_LONG, self_rank+stride, MPI_ANY_TAG, communicator, &recv_req);
            MPI_Wait(&recv_req , &recv_status);
            local_sum += recv_buf; // perform pairwise sum here
            printf("rank: %d received\n", self_rank);
        }
        else{ // odd ranks after stride: sender
            MPI_Isend(&local_sum, 1, MPI_LONG_LONG, self_rank-stride, 0, communicator, &send_req);
            MPI_Wait(&send_req , &send_status);
            printf("rank: %d sent\n", self_rank);
        }
        stride *= 2;
        // printf("rank: %d about to hit MPI_Barrierr\n\n", self_rank);
        MPI_Barrier(communicator); // sync here
        // printf("rank: %d just hit barrier\n\n", self_rank);
    }
    if (self_rank == root){*recv_data = local_sum;}
    return 0;
}




int main(int argc, char* argv[]){
  // Initialize the MPI environment
    int world_rank, world_size; // init world rank and size
    MPI_Init(&argc, &argv);

    // Find out rank, size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // size of a array is determined by how many nodes are working on this task
    int arrSize = 2<<30 / world_rank; 

    long long int* bigArr = malloc(sizeof(long long int)*arrSize);

    for (int i=0; i<arrSize; i++){
        bigArr[i] = world_rank * arrSize;
    }

    // calling MPI_P2P_Reduce
    long long int global_sum = 0; // only 0th rank process will have this value modified 
    uint64_t p2p_start_cycles = clock_now();
    MPI_P2P_Reduce(bigArr, &global_sum, arrSize, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    uint64_t p2p_end_cycles = clock_now();

    if (world_rank == 0){ // only root has the result and it prints it here
        long long int correct_answer = 576460751766552576;
        if (correct_answer - global_sum == 0){
            printf("result is correct: %lld\n", correct_answer);
        }
        else{
            printf("correct answer: %lld\nmy answer:      %lld\n",correct_answer, global_sum);
        }
    }

    //calling MPI_Reduce


    // show runtime
    if (world_rank == 0){
        double p2p_time_in_secs = ((double)(p2p_end_cycles - p2p_start_cycles)) / 512000000;
        printf("MPI_P2P_Reduce took %f seconds.\n", p2p_time_in_secs);
    }

    free(bigArr);
    MPI_Finalize();
    return 0;
}
