#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "clockcycle.h"


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
int main(int argc, char *argv[]) {
  
  int npes, myrank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  long long arrSize = (1 << 30) / npes;
  long long * bigArr = malloc(arrSize * sizeof(long long));
 
  for (int i = 0; i < arrSize; i++) { 
    bigArr[i] = i + (long long) myrank * arrSize;
  }
  unsigned long long local_sum = 0;
  for (int i = 0; i < arrSize; i++) { 
      local_sum += bigArr[i];
  }
  long long rootSum;
  unsigned long long start_time, end_time, reduced_sum;
  start_time = clock_now();
  MPI_P2P_Reduce(bigArr, &rootSum, arrSize, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  end_time = clock_now();


  unsigned long long start_time2, end_time2, reduced_res;
  start_time2 = clock_now();
  MPI_Reduce(&local_sum, &reduced_sum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  end_time2 = clock_now();
  
  
  if (myrank == 0) {
    printf("Reduce Completed in %f seconds. MPI_Reduce completed in %f seconds\n",
            (double) (end_time - start_time) / 512e6, (double) (end_time2 - start_time2) / 512e6);
  }
  free(bigArr);
  MPI_Finalize();
  if (myrank == 0) {
      printf("p2p reduce: %lld MPI_Reduce: %lld\n", rootSum, reduced_sum);
  }
  return EXIT_SUCCESS;
}