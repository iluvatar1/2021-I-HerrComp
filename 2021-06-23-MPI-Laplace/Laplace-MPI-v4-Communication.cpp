#include <iostream>
#include <vector>
#include <mpi.h>

const int NX = 12;
const int NY = 12;
const double XMIN = 0.0;
const double XMAX = 1.0;
const double YMIN = 0.0;
const double YMAX = 1.0;
const double DELTAX = (XMAX-XMIN)/NX;
const double DELTAY = (YMAX-YMIN)/NY;
const int NSTEPS = 100;

typedef std::vector<double> data_t;

// serial functions
void initial_conditions(data_t & data, int nx, int ny);
void boundary_conditions(data_t & data, int nx, int ny);
void evolve(data_t & data, int nx, int ny, int nsteps);
void relaxation_step(data_t & data, int nx, int ny);
void print_screen(const data_t & data, int nx, int ny);
void start_gnuplot(void);
void print_gnuplot(const data_t & data, int nx, int ny);

// parallel functions
void mpi_print_screen(const data_t & data, int nx, int ny,
                      int pid, int np);
void mpi_print_screen_aux(const data_t & data, int nx, int ny);
void mpi_initial_conditions(data_t & data, int nx, int ny,
                            int pid, int np);
void mpi_communicate(data_t & data, int nx, int ny,
                     int pid, int np);

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int pid = 0;
    int np = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    // declare data structures
    int Nl = NX/np + 2;
    data_t potential(Nl*NY); // [ii, jj] -> ii*NY + jj
    //mpi_print_screen(potential, Nl, NY, pid, np);

    // set initial and boundary conditions
    //initial_conditions(potential, NX, NY);
    // boundary_conditions(potential, NX, NY);
    mpi_initial_conditions(potential, Nl, NY, pid, np);
    mpi_print_screen(potential, Nl, NY, pid, np);
    mpi_communicate(potential, Nl, NY, pid, np);
    mpi_print_screen(potential, Nl, NY, pid, np);

    // evolve and print
    //evolve(potential, NX, NY, NSTEPS);

    MPI_Finalize();

    return 0;
}

// parallel functions
void mpi_print_screen(const data_t & data, int nx, int ny,
                      int pid, int np)
{
    if (0 == pid) {
        data_t buffer(nx*ny);
        mpi_print_screen_aux(data, nx, ny);
        // recibo de los demas e imprimo
        for (int id = 1; id < np; ++id) {
            MPI_Recv(&buffer[0], nx*ny, MPI_DOUBLE, id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mpi_print_screen_aux(buffer, nx, ny);
        }
    } else {
        MPI_Send(&data[0], nx*ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

void mpi_print_screen_aux(const data_t & data, int nx, int ny)
{
    for(int ix = 0; ix < nx; ++ix) {
        if (0 == ix || nx-1 == ix) std::cout << "G: ";
        else std::cout << "   ";
        for(int iy = 0; iy < ny; ++iy) {
            std::cout << data[ix*ny + iy] << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}
void mpi_initial_conditions(data_t & data, int nx, int ny, int pid, int np)
{
    for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
            data[ix*ny + iy] = pid;
        }
    }
}

void mpi_communicate(data_t & data, int nx, int ny,
                     int pid, int np)
{
    // forward send
    if (pid != np-1) {
        MPI_Send(&data[(nx-2)*ny], ny, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD);
    }
    if (0 != pid) {
        MPI_Recv(&data[0], ny, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

// Serial functions
void initial_conditions(data_t & data, int nx, int ny)
{
    for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
            data[ix*ny + iy] = 1.0;
        }
    }
}
void boundary_conditions(data_t & data, int nx, int ny)
{
    int ix, iy;
    // first row
    ix = 0;
    for(int iy = 0; iy < ny; ++iy) {
        data[ix*ny + iy] = 100.0;
    }
    // last row
    ix = nx-1;
    for(int iy = 0; iy < ny; ++iy) {
        data[ix*ny + iy] = 0.0;
    }
    // first column
    iy = 0;
    for(int ix = 1; ix < nx; ++ix) {
        data[ix*ny + iy] = 0.0;
    }
    // last column
    iy = ny-1;
    for(int ix = 1; ix < nx; ++ix) {
        data[ix*ny + iy] = 0.0;
    }
    //new
    //ix = nx/2;
    //for(int iy = ny/3; iy <= 2*ny/3; ++iy) {
    //    data[ix*ny + iy] = -50.0;
    //}
}

void evolve(data_t & data, int nx, int ny, int nsteps)
{
    start_gnuplot();
    for(int istep = 0; istep < nsteps; ++istep) {
        relaxation_step(data, nx, ny);
        //print_screen(data, nx, ny);
        print_gnuplot(data, nx, ny);
    }
}
void relaxation_step(data_t & data, int nx, int ny)
{
    // recorrer toda la matriz y aplicar el algoritmo,
    // teniendo cuidado con no modificar las condiciones de
    // frontera
    for(int ix = 1; ix < nx-1; ++ix) {
        for(int iy = 1; iy < ny-1; ++iy) {
            // check that this cell is NOT a boundary condition or a border
            //if ( (ix == nx/2) && (ny/3 <= iy) && (iy <= 2*ny/3) ) continue;
            // update the cell
            data[ix*ny + iy] = (data[(ix+1)*ny + iy] + data[(ix-1)*ny + iy] + data[ix*ny + iy+1] + data[ix*ny + iy-1])/4.0;
        }
    }

}

void print_screen(const data_t & data, int nx, int ny)
{
    for(int ix = 0; ix < nx; ++ix) {
        for(int iy = 0; iy < ny; ++iy) {
            std::cout << data[ix*ny + iy] << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void start_gnuplot(void)
{
    std::cout << "set pm3d\n";
    std::cout << "set contour base\n";
    std::cout << "set term gif animate\n";
    std::cout << "set output 'anim.gif'\n";
}

void print_gnuplot(const data_t & data, int nx, int ny)
{
    std::cout << "splot '-' w l lt 3 \n";
    for(int ix = 0; ix < nx; ++ix) {
        double x = XMIN + ix*DELTAX;
        for(int iy = 0; iy < ny; ++iy) {
            double y = YMIN + iy*DELTAY;
            std::cout << x << "  " << y << "  " << data[ix*ny + iy] << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "e\n";
}
