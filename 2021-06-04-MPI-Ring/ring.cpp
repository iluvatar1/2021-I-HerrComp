#include <iostream>
#include <cstdlib>
#include <cmath>
#include "mpi.h"

void ring(int nreps, int pid, int np);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv); /* Mandatory */

  int pid;                 /* rank of process */
  int np;                 /* number of processes */

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  int NREPS = std::atoi(argv[1]);

  ring(NREPS, pid, np);

  MPI_Finalize(); /* Mandatory */
  
  return 0; 
}

void ring(int nreps, int pid, int np)
{
  int tag = 0;
  double totaltime = 0;
  int val = pid;
  int buf = pid;
  for (int isteps = 0; isteps < nreps; ++isteps){
    val = pid;
    buf = pid;
    int next = (pid+1)%np;
    int prev = (pid-1+np)%np;
    double start = MPI_Wtime();
    if (pid == 0) {
      MPI_Send(&pid, 1, MPI_INT, next, tag, MPI_COMM_WORLD);
      MPI_Recv(&buf, 1, MPI_INT, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else { // pid != 0
      MPI_Recv(&buf, 1, MPI_INT, prev, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      val += buf;
      MPI_Send(&val, 1, MPI_INT, next, tag, MPI_COMM_WORLD);
    }
    double end = MPI_Wtime();
    totaltime += end-start;
  }
  if (pid == 0) {
    //std::cout << "totaltime: " << totaltime  << std::endl;
    //std::cout << "totaltime/nsteps: " << totaltime/nsteps  << std::endl;
    std::cout << "bandwidth: " << 2.0*nreps*1*sizeof(double)/totaltime/1.0e9  << std::endl;
    std::cout << "val: " << buf << std::endl;
  }
}
