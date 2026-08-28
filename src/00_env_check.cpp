#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fftw3.h>

#include <cmath>
#include <complex>
#include <iostream>

static void check_eigen(){
    const int n = 5;
    Eigen::SparseMatrix<double> A(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i=0; i<n; ++i){
        triplets.emplace_back(i, i, 2.0);
        if (i > 0) triplets.emplace_back(i, i-1, -1.0);
        if (i < n-1) triplets.emplace_back(i, i+1, -1.0);
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd x = solver.solve(b);

    std::cout << "[Eigen] tridiagonal solve residual = "
              << (A*x - b).norm() << "\n";
}

int main(){
    std::cout << "Env check: CMake + Eigen\n";
    check_eigen();
    return 0;

}