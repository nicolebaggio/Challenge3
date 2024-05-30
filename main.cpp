#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <omp.h>



/*!
    * \brief Function to calculate the value of the function f(x, y) = 8 * pi^2 * sin(2 * pi * x) * cos(2 * pi * y)
    * \param x The x coordinate
    * \param y The y coordinate
    * \return The value of the function at (x, y)
    
*/
const double pi = 3.14159265358979323846;
std::function<double(double, double)> f = [](double x, double y) -> double {
    return 8 * pi * pi * sin(2 * pi * x) * cos(2 * pi * y);
};

/*!
    * \brief Function to check the convergence of the solution
    * \param U The current solution
    * \param U1 The previous solution
    * \param h The step size
    * \return 1 if the solution has converged, 0 otherwise
*/
int check_convergence(const std::vector<std::vector<double>>& U, const std::vector<std::vector<double>>& U1, double h) {
    double tol = 1e-1;
    double s = 0.0;
    #pragma omp parallel for reduction(+:s)
    for (size_t i = 0; i < U.size(); i++) {
        for (size_t j = 0; j < U[i].size(); j++) {
            s += (U[i][j] - U1[i][j]) * (U[i][j] - U1[i][j]);
        }
    }
    if (std::sqrt(h * s) < tol) {
        return 1;
    }
    return 0;
};


/*!
    * \brief Function to write the solution to a VTK file
    * \param U The solution
    * \param n The size of the matrix
    * \param h The step size
    * \param filename The name of the file
*/
void write_vtk(const std::vector<std::vector<double>>& U, int n, double h, const std::string& filename) {
    std::ofstream file(filename);
    file << "# vtk DataFile Version 2.0\n";
    file << "Solution\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << n << " " << n << " 1\n";
    file << "ORIGIN 0 0 0\n";
    file << "SPACING " << h << " " << h << " " << h << "\n";
    file << "POINT_DATA " << n * n << "\n";
    file << "SCALARS solution double 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file << U[i][j] << "\n";
        }
    }

    file.close();
}



/*!
    * \brief Main function
    * \param argc The number of command line arguments
    * \param argv The command line arguments
    * \return 0 if the program runs successfully, 1 otherwise
*/
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int n,size;

    if (rank == 0) {
        std::cout << "Enter the size of the matrix:" << std::endl;
        std::cin >> n;
        if (n < 1) {
            std::cout << "The size of the matrix should be greater than 0" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        std::cout << "Enter the number of processor you want to use:" << std::endl;
        std::cin >> size;
        if (size < 1 || size > world_size) {
            std::cout << "The number of processors should be between 1 and " << world_size << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        for (int i = 1; i < world_size; i++) {
            MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&size, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

    }
    else {
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank >= size) {
        MPI_Finalize();
        return 0;
    }


    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U1(n, std::vector<double>(n, 0.0));

    std::vector<int> local_r(size, 0);
    std::vector<int> local_start_idx(size, 0);

    int start_idx = 0;
    for (int i = 0; i < size; i++) {
        local_r[i] = (n / size) + (i < (n % size) ? 1 : 0);
        local_start_idx[i] = start_idx;
        start_idx += local_r[i];
    }

    double h = 1.0 / (n - 1);

    for (int iter = 0; iter < 100; iter++) {
        // Synchronize the boundary rows between processes
        if (rank > 0) {
            MPI_Recv(U[local_start_idx[rank] - 1].data(), n, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Send(U[local_start_idx[rank] + local_r[rank] - 1].data(), n, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }

        // Update the matrix in each process
        #pragma omp parallel for collapse(2)
        for (int i = local_start_idx[rank]; i < local_start_idx[rank] + local_r[rank]; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                if (i > 0 && i < n - 1) {
                    U1[i][j] = 0.25 * (U[i-1][j] + U[i+1][j] + U[i][j-1] + U[i][j+1] + h * h * f(i * h, j * h));
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // Check for convergence
        int local_convergence = check_convergence(U, U1, h);
        int global_convergence;
        MPI_Allreduce(&local_convergence, &global_convergence, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

        // Update U with U1
        #pragma omp parallel for
        for (int i = local_start_idx[rank]; i < local_start_idx[rank] + local_r[rank]; ++i) {
            for (int j = 0; j < n; ++j) {
                U[i][j] = U1[i][j];
            }
        }

        if (global_convergence == 1) {
            break;
        }
    }

    // Gather the full matrix on rank 0
    if (rank == 0) {
        for (int r = 1; r < size; ++r) {
            for (int i = local_start_idx[r]; i < local_start_idx[r] + local_r[r]; ++i) {
                MPI_Recv(U[i].data(), n, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    } else {
        for (int i = local_start_idx[rank]; i < local_start_idx[rank] + local_r[rank]; ++i) {
            MPI_Send(U[i].data(), n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    // Print the matrix on rank 0
    if (rank == 0) {
       // for (int i = 0; i < n; i++) {
           // for (int j = 0; j < n; j++) {
             //   std::cout << U[i][j] << " ";
         //   }
          //  std::cout << std::endl;
     //   }
        write_vtk(U, n, h, "solution.vtk");
        std::cout << "Solution written to solution.vtk" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
