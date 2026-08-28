#include "fd/grid1d.hpp"
#include "fd/poisson1d/assemble.hpp"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
#include <vector>
using namespace fd;
using namespace fd::poisson1d;

int main(){
    const double pi = std::acos(-1.0); 
    const double a = 0.0;
    const double b = 2*pi;
    std::vector<int> Ns = {10, 50, 100, 500, 1000};
    for (const int N : Ns){
        Grid1D grid(a, b, N);
        Eigen::SparseMatrix<double> A = assembleA(grid);
        Eigen::MatrixXd dense(A);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(dense);
        double eig = solver.eigenvalues()(0); //smallest eigenvalue
        std::cout<< "N = " << grid.N() <<"; h = "<< grid.h() <<"; smallest eigenvalue = " << eig <<"\n";
    }
    std::cout<<"Eigenvalue Limit = " << (pi/(b-a))*(pi/(b-a))<<"\n";
}