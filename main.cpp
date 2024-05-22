#include <vector>
#include <mpi.h>
#include <omp.h>
#include <functional>
#include <iostream>
#include <cmath>


const double pi=3.14159265358979323846;
std::function<double(double , double )> f=[](double x,double y)->double{
    return 8*pi*pi* sin(pi*x)*cos(pi*y);
};


int check_convergence(std::vector<std::vector<double>> U, std::vector<std::vector<double>> U1,double h);

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   

    int n;
    int size;
    std::vector<std::vector<double>> U(0);

    if(rank==0){
        std::cout<<"Enter the number of parallel tasks"<<std::endl;
        std::cin>>size;
        std::cout<<"Enter the size of the matrix"<<std::endl;
        std::cin>>n;
    }
    

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    U.resize(n);
    for(int i=0; i<n; i++){
        U[i].resize(n,0.0);
    }; 
   

    int local_nrows = (n % size > rank) ? n/size +1 : n/size;

    int start_idx = 0;
    std::vector<int> recv_counts(size);
    std::vector<int> send_counts(size);
    std::vector<int> recv_start_idx(size);
    std::vector<int> send_start_idx(size);

  for (int i = 0; i < size; ++i)
    {
      recv_counts[i] = (n % size > i) ? n / size + 1 : n / size;
      send_counts[i] = recv_counts[i] * n;

      recv_start_idx[i] = start_idx;
      send_start_idx[i] = start_idx * n;

      start_idx += recv_counts[i];
    }

    int h=1/n;
    std::vector<std::vector<double>> U1=U;

    int conv_tot=0;

    if(rank==0){
      for(int iter=0; iter<100; iter++){
      for (int i=1; i<recv_start_idx[i]; ++i){
        for(int j=1; j<n; ++j){
          U1[i][j]=1/4 * (U[i-1][j]+U[i+1][j]+U[i][j-1]+U[i][j+1]+h*h*f(i*h,j*h));
        }
      }
      MPI_Send(U[recv_start_idx[0]+recv_counts[0]-1].data(), n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
      int convergence=check_convergence(U,U1,h);
      int conv_tot=0;
      MPI_Allgather(&convergence, 1, MPI_INT, &conv_tot, 1, MPI_INT, MPI_COMM_WORLD);
      if(conv_tot==size){
            break;
      }
      for (int i=1; i<recv_start_idx[i]; ++i){
        for(int j=1; j<n; ++j){
          U[i][j]=U1[i][j];
        }
        MPI_Gather(U[i].data(), n, MPI_DOUBLE, U[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      } 
    }

    }


    for (int i=1; i<size; ++i){
      if(rank==i){
        for(int iter=0; iter<100; iter++){
        std::vector<double> previous_vec(n,0.0);
        MPI_Recv(U[recv_start_idx[i]-1].data(), n, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i=recv_start_idx[i]; i<recv_start_idx[i]+recv_counts[i]; ++i){
          for(int j=1; j<n; ++j){
            U1[i][j]=1/4 * (U[i-1][j]+U[i+1][j]+U[i][j-1]+U[i][j+1]+h*h*f(i*h,j*h));
          }
        }
        if(i!=size)
            MPI_Send(U[recv_start_idx[i]+recv_counts[i]-1].data(), n, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD);
        int convergence=check_convergence(U,U1,h);
        int conv_tot=0;
        MPI_Allgather(&convergence, 1, MPI_INT, &conv_tot, 1, MPI_INT, MPI_COMM_WORLD);
        if(conv_tot==size){
            break;
        }
        for (int i=1; i<recv_start_idx[i]; ++i){
           for(int j=1; j<n; ++j){
              U[i][j]=U1[i][j];
            }
        MPI_Gather(U[recv_start_idx[i]].data(), n, MPI_DOUBLE, U[recv_start_idx[i]].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }  

      }
      }
    }

  
    
    MPI_Finalize();
    return 0;
}


int check_convergence(std::vector<std::vector<double>> U, std::vector<std::vector<double>> U1, double h){
    double tol=1e-1;
    double s=0.0;
    for(int i=0; i<U.size(); i++){
        for(int j=0; j<U[i].size(); j++){
               s+=abs((U[i][j]-U1[i][j])*(U[i][j]-U1[i][j]));
        }
    }
    if(std::sqrt(h*s)<tol){
        return 1;
    }
    return 0;
}