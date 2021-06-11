#include <iostream>
#include <cstdlib>
#include <cmath>
#include "mpi.h"

void average(int nsize, int pid, int np);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv); /* Mandatory */

  int pid;                 /* rank of process */
  int np;                 /* number of processes */

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  int NSIZE = std::atoi(argv[1]);
  double start = 0; 

  /*
  if (0 == pid) {
    std::cout << "Escribe el tamanho del arreglo:\n";
    std::cin >> NSIZE;
    start = MPI_Wtime();
    for(int npid = 1; npid < np; ++npid) {
      MPI_Send(&NSIZE, 1, MPI_INT, npid, 0, MPI_COMM_WORLD);
    }
  } else {
      MPI_Recv(&NSIZE, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
  }
  */

  /*if (0 == pid) {
    std::cout << "Escribe el tamanho del arreglo:\n";
    std::cin >> NSIZE;
    }*/
  start = MPI_Wtime();
  MPI_Bcast(&NSIZE, 1, MPI_INT, 0, MPI_COMM_WORLD);
  double end = MPI_Wtime();

  if (0 == pid) {
    std::cout << "Comm time: " << end-start << "\n";
  }
  
  average(NSIZE, pid, np);

  MPI_Finalize(); /* Mandatory */
  
  return 0; 
}

void average(int nsize, int pid, int np)
{
  int tag1 = 0, tag2 = 1;

  int Nlocal = nsize/np;
  double *data = new double [Nlocal] {0};

  // falta inicializar el arreglo!
  for (int ii = 0; ii < Nlocal; ++ii) {
    data[ii] = 1;
  }
  
  double suma = 0;
  for (int ii = 0; ii < Nlocal; ++ii) {
    suma += data[ii];
  }

  // communicate results
  double total = 0;
  MPI_Reduce(&suma, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (0 == pid) {
    std::cout << "El promedio es: " << total/nsize << "\n";
  }
  /*
  if (0 == pid) {
    double tmp = 0;
    for(int npid = 1; npid < np; ++npid) {
      MPI_Recv(&tmp, 1, MPI_DOUBLE, npid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      suma += tmp;
    }
    std::cout << "El promedio es: " << suma/nsize << "\n";
  } else {
      MPI_Send(&suma, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
  */
  delete [] data;
}
