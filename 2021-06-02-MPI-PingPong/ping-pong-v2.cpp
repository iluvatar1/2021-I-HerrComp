#include <iostream>
#include <cstdlib>
#include <cmath>
#include "mpi.h"

void pingpong(int nsize, int nsteps, int pid, int np);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv); /* Mandatory */

  int pid;                 /* rank of process */
  int np;                 /* number of processes */

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  int NSIZE = std::atoi(argv[1]);
  int NSTEPS = std::atoi(argv[2]);

  pingpong(NSIZE, NSTEPS, pid, np);

  MPI_Finalize(); /* Mandatory */
  
  return 0; 
}

void pingpong(int nsize, int nsteps, int pid, int np)
{
  int tag1 = 0, tag2 = 1;

  double totaltime = 0;
  double *data = new double [nsize] {0};
  for (int isteps = 0; isteps < nsteps; ++isteps){
    double start = MPI_Wtime();
    if (pid == 0) {
      MPI_Send(data, nsize, MPI_DOUBLE, 1, tag1, MPI_COMM_WORLD);
      MPI_Recv(data, nsize, MPI_DOUBLE, 1, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else { // pid == 1
      MPI_Recv(data, nsize, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(data, nsize, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD);
    }
    double end = MPI_Wtime();
    totaltime += end-start;
  }
  if (pid == 0) {
    //std::cout << "totaltime: " << totaltime  << std::endl;
    //std::cout << "totaltime/nsteps: " << totaltime/nsteps  << std::endl;
    std::cout << "bandwidth: " << 2.0*nsteps*nsize*sizeof(double)/totaltime/1.0e9  << std::endl;
  }
  delete [] data;
}
