#include <iostream>
#include <cstdlib>
#include <cmath>
#include "mpi.h"

void scatter(int nsize, int pid, int np);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv); /* Mandatory */

  int pid;                 /* rank of process */
  int np;                 /* number of processes */

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  int NSIZE = std::atoi(argv[1]);

  scatter(NSIZE, pid, np);

  MPI_Finalize(); /* Mandatory */
  
  return 0; 
}

void scatter(int nsize, int pid, int np)
{
  double *data = nullptr;
  int nlocal = nsize/np;
  double *buf = new double [nlocal] {0};

  if (0 == pid) {
    data = new double [nsize] {0};
    // falta inicializar el arreglo!
    for (int ii = 0; ii < nsize; ++ii) {
      data[ii] = (2*ii+1);
      //std::cout << data[ii] << std::endl;
    }
  }
  MPI_Scatter(data, nlocal, MPI_DOUBLE, buf, nlocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  std::cout << "pid: " << pid << " -> ";
  for (int ii = 0; ii < nlocal; ++ii) {
    std::cout << buf[ii] << "  ";
  }
  std::cout << "\n";
  
  delete [] data;
  delete [] buf;
}
