#include <iostream>
#include <cstdlib>
#include "mpi.h"

using fptr = double(double);

double f(double x);
double integral(double a, double b, int nsteps, fptr fun, int pid, int np);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv); /* Mandatory */

  int pid;                 /* rank of process */
  int np;                 /* number of processes */

  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  double a = std::atof(argv[1]);
  double b = std::atof(argv[2]);
  int N = std::atoi(argv[3]);;
  
  std::cout << integral(a, b, N, f, pid, np) << std::endl;
  
  MPI_Finalize(); /* Mandatory */
  
  return 0;
}

double f(double x)
{
  return x;
}

double integral(double a, double b, int nsteps, fptr fun, int pid, int np)
{
  double result = 0;
  double dx = (b-a)/nsteps;
  int Nlocal = nsteps/np;
  int iimin = pid*Nlocal;
  int iimax = (pid+1)*Nlocal;
  for (int ii = iimin; ii < iimax; ++ii) {
    double x = a + ii*dx;
    result += fun(x);
  }
  result *= dx;
  return result;
}
