#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "clockcycle.h"


// this function takes the same inputs as MPI_reduce, except MPI_Op is set to MPI_SUM
void MPI_P2P_Reduce(long long *sendbuf, void *recvbuf, int count,
                    MPI_Datatype datatype, int target,
                    MPI_Comm comm) {
  int stride = 1;
  int rank, size;
  int curRank = MPI_Comm_rank(comm, &rank);

  long long sum = 0; 
  for (int i = 0; i < count; i++) { 
    sum += sendbuf[i];
  }
  MPI_Barrier(comm);
  //printf("count: %d sum: %lld\n", count, sum);
  //printf("%d", MPI_Comm_size(comm, &size));
  MPI_Comm_size(comm, &size);
  //printf("size: %d\n", size);
  //MPI_Comm_size(comm, &size);
  while (stride < size) {
    //printf("%d %d\n", size, rank);       
    if ((rank / stride) % 2) {  // sender
      MPI_Request req;
      MPI_Isend(&sum, 1, MPI_LONG_LONG, rank - stride, 0, comm, &req);
      MPI_Wait(&req, MPI_STATUS_IGNORE);
    } else {  // receiver
      long long temp;
      MPI_Request req;
      MPI_Irecv(&temp, 1, MPI_LONG_LONG, rank + stride, MPI_ANY_TAG, comm, &req);
      MPI_Wait(&req, MPI_STATUS_IGNORE);
      sum += temp;
      if (sum >= 576460751766552576) 
          printf("temp: %lld\n", sum);
    }
    stride *= 2;  
    MPI_Barrier(comm);
  }
  return;
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
