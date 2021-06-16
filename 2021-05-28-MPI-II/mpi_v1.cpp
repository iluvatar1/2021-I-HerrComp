#include <iostream>
#include "mpi.h"

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv); /* Mandatory */

  int pid;                 /* rank of process */
  int np;                 /* number of processes */

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  for (int id = 0; id < np; id++) {
    if (id == pid) 
      std::cout << "Hello from pid: " << pid << " out of np: " << np << "\n";
  }
  
  MPI_Finalize(); /* Mandatory */
  
  return 0;
}
