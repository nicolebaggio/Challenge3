#include <vector>
#include <mpi.h>
#include <omp.h>
#include <functional>
#include <iostream>
#include <cmath>


const double pi=3.14159265358979323846;
std::function<double(double , double )> f=[](double x,double y)->double{
    return 8*pi*pi* sin(2*pi*x)*cos(2*pi*y);
};


int check_convergence(std::vector<std::vector<double>> U, std::vector<std::vector<double>> U1,double h){
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
};

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   

    int n,size;

    if(rank==0){
        std::cout<<"Enter the number of parallel tasks:"<<std::endl;
        std::cin>>size;
        if (size<1){
            std::cout<<"The number of parallel tasks should be greater than 0"<<std::endl;
            return 0;
        }
        std::cout<<"Enter the size of the matrix:"<<std::endl;
        std::cin>>n;
        if (n<1){
            std::cout<<"The size of the matrix should be greater than 0"<<std::endl;
            return 0;
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<std::vector<double>> U(n,std::vector<double>(n,0.0));
    std::vector<std::vector<double>> U1(n,std::vector<double>(n,0.0));

    std::vector<int> local_r(size,0);
    std::vector<int> local_start_idx(size,0);

    int start_idx=0;
    for(int i=0; i<size; i++){
        local_r[i]=(n % size<i) ? n/size+1 : n/size;
        local_start_idx[i]=start_idx;
        start_idx+=local_r[i];
    }

    double h=1.0/(n-1);

    if(rank==0){
      for(int iter=0; iter<100; iter++){
        for (int i=1; i<local_r[0]; ++i){
          for(int j=1; j<n; ++j){
            U1[i][j]=1/4 * (U[i-1][j]+U[i+1][j]+U[i][j-1]+U[i][j+1]+h*h*f(i*h,j*h));
          }
      }

      MPI_Send(U1[local_r[0]-1].data(), n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
      int convergence=check_convergence(U,U1,h);
      for (int i=local_start_idx[0]; i<local_start_idx[1]; ++i){
             for(int j=1; j<n; ++j){
               U[i][j]=U1[i][j];
              }
          }  
       if(convergence==1){
         for(int i=local_start_idx[r]; i<local_start_idx[r+1]; ++i){
            MPI_Gather(U[i].data(), n, MPI_DOUBLE, U[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         }
           break;
      }
    }
    }

    for (int r=1; r<size; ++r){
      if(rank==r && r!=size-1){
        for(int iter=0; iter<100; iter++){
           MPI_Barrier(MPI_COMM_WORLD); 
           MPI_Recv(U[local_start_idx[r]-1].data(), n, MPI_DOUBLE, r-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
           for (int i=local_start_idx[r]; i<local_start_idx[r+1]; ++i){
             for(int j=1; j<n; ++j){
               U1[i][j]=1/4 * (U[i-1][j]+U[i+1][j]+U[i][j-1]+U[i][j+1]+h*h*f(i*h,j*h));
             }
           }
        
           MPI_Send(U1[local_start_idx[r]+local_r[r]-1].data(), n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
    
          int convergence=check_convergence(U,U1,h);
          for (int i=local_start_idx[r]; i<local_start_idx[r+1]; ++i){
             for(int j=1; j<n; ++j){
               U[i][j]=U1[i][j];
              }
          }  
          if(convergence==1){
            for(int i=local_start_idx[r]; i<local_start_idx[r+1]; ++i){
              MPI_Gather(U[i].data(), n, MPI_DOUBLE, U[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            break;
         }

      }
      }
      else if(rank==r && r==size-1){
        for(int iter=0; iter<100; iter++){
           MPI_Recv(U[local_start_idx[r]-1].data(), n, MPI_DOUBLE, r-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
           for (int i=local_start_idx[r]; i<local_start_idx[r]+local_r[r]; ++i){
             for(int j=1; j<n; ++j){
               U1[i][j]=1/4 * (U[i-1][j]+U[i+1][j]+U[i][j-1]+U[i][j+1]+h*h*f(i*h,j*h));
             }
           }
           MPI_Barrier(MPI_COMM_WORLD);    
           int convergence=check_convergence(U,U1,h);
           for (int i=local_start_idx[r]; i<local_start_idx[r+1]; ++i){
             for(int j=1; j<n; ++j){
               U[i][j]=U1[i][j];
              }
          }  
          if(convergence==1){
            for(int i=local_start_idx[r]; i<local_start_idx[r+1]; ++i){
              MPI_Gather(U[i].data(), n, MPI_DOUBLE, U[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            break;
         }
        }
      }

    }
   
    


    if (rank == 0) {
       for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cout << U[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    
    MPI_Finalize();
    return 0;
}

