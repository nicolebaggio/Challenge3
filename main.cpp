#include <vector>
#include <mpi.h>
#include <omp.h>
#include <functional>
#include <iostream>

std::function<double(double , double )> f=[](double x,double y)->double{
    return 8*pi*pi* sin(pi*x)*cos(pi*y);
}

int main(*argc, *argv) {
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
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    U.resize(n);
    for(int i=0; i<n; i++){
        U[i].resize(n,0.0);
    }; 
   

    int local_nrows = (n % size > rank) ? n/size +1 : n/size;

    int start_idx = 0;
  for (int i = 0; i < size; ++i)
    {
      recv_counts[i] = (nrows % size > i) ? nrows / size + 1 : nrows / size;
      send_counts[i] = recv_counts[i] * ncols;

      recv_start_idx[i] = start_idx;
      send_start_idx[i] = start_idx * ncols;

      start_idx += recv_counts[i];
    }

    int h=1/n;
    std::vector<std::vector<double>> U1=U;

    int conv_tot=0;

    if(rank==0){
      for (int i=1; i<recv_start_idx[i]; ++i){
        for(int j=1; j<n; ++j){
          U1[i][j]=1/4 * (U[i-1][j]+U[i+1][j]+U[i][j-1]+U[i][j+1]+h*h*f(i,j));
        }
      }
      MPI_Send(U[recv_start_idx[0]+recv_counts[0]-1].data(), n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
      bool convergence=check_convergence(U,U1,h);
        MPI_All_Gather(&convergence, 1, MPI_INT, conv_tot, 1, MPI_INT, MPI_COMM_WORLD);
        if(conv_tot==size){
            break;
        }
      for (int i=1; i<recv_start_idx[i]; ++i){
        for(int j=1; j<n; ++j){
          U[i][j]=U1[i][j];
        }
      }  
    }


    for (int i=1; i<size; ++i){
      if(rank==i){
        std::vector<double> previous_vec(n,0.0);
        MPI_Recv(U[recv_start_idx[i]-1].data(), n, MPI_DOUBLE, i-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i=recv_start_idx[i]; i<recv_start_idx[i]+recv_counts[i]; ++i){
          for(int j=1; j<n; ++j){
            U1[i][j]=1/4 * (U[i-1][j]+U[i+1][j]+U[i][j-1]+U[i][j+1]+h*h*f(i,j));
          }
        }
        if(i!=size)
            MPI_Send(U[recv_start_idx[i]+recv_counts[i]-1].data(), n, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD);
        bool convergence=check_convergence(U,U1,h);
        MPI_All_Gather(&convergence, 1, MPI_INT, conv_tot, 1, MPI_INT, MPI_COMM_WORLD);
        if(conv_tot==size){
            break;
        }
        for (int i=1; i<recv_start_idx[i]; ++i){
        for(int j=1; j<n; ++j){
          U[i][j]=U1[i][j];
        }
      }  

      }
    }

  
    
    MPI_Finalize();
    return 0;
}


bool check_convergence(std::vector<std::vector<double>> U, std::vector<std::vector<double>> U1,double h){
    double tol=1e-1;
    double s=0.0;
    for(int i=0; i<U.size(); i++){
        for(int j=0; j<U[i].size(); j++){
            if(abs(U[i][j]-U1[i][j])>max){
               s+=abs((U[i][j]-U1[i][j])*(U[i][j]-U1[i][j]));
            }
        }
    }
    if(std::sqrt(h*s)<tol){
        return true;
    }
    return false;
}