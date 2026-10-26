#include "fd/grid1d.hpp"
#include "fd/poisson1d/boundary_condition.hpp"
#include "fd/sparse_solve.hpp"
#include "core/manufactured_solutions.hpp"
#include "fd/poisson1d/scheme.hpp"
#include "core/convergence.hpp"
#include "fd/error_norms.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>
#include <vector>

using namespace core;
using namespace fd;
using namespace fd::poisson1d;

int main(){
    /* Setup 1d Poisson for a test case and
    compare against manufactured solutions*/
    std::vector<int> Ns = {5, 10, 20, 50, 100, 200, 500, 1000, 2000};
    int n = Ns.size();
    std::function<double(double)> rho_fn = [](double x){return 2*std::sin(x)+x*std::cos(x);};
    std::function<double(double)> phi_fn = [](double x){return x*std::cos(x);};
    ManufacturedSolution ms{rho_fn, phi_fn};
    // steup domain and mesh
    const double a = 0.;
    const double b = 2*std::acos(-1.0); //2pi
    // bc values
    const double alpha = phi_fn(a);
    const double beta = phi_fn(b);
    DirichletBC bc(alpha, beta);
    const SecondOrderScheme scheme;
    Eigen::VectorXd H(n);
    Eigen::VectorXd l1_errs(n);
    Eigen::VectorXd l2_errs(n);
    Eigen::VectorXd linf_errs(n);
    std::ofstream file("data/poisson_convergence.csv");
    file << "N,h,l1,l2,linf\n";
    int j = 0;
    for (const int N : Ns){
        const Grid1D grid(a, b, N);
        const double h = grid.h();
        Eigen::SparseMatrix<double> A = bc.assembleA(grid, scheme);
        Eigen::VectorXd rhs = bc.assembleRHS(grid, scheme, ms.rho);
        Eigen::VectorXd reduced_sol = solveSPD(A, rhs);
        Eigen::VectorXd sol = bc.reconstructSolution(grid, reduced_sol);
        // exact solution
        Eigen::VectorXd sol_exact(grid.N()+1);
        Eigen::VectorXd nodes = grid.nodes();
        for (int i = 0; i<N+1; i++){
            sol_exact(i) = ms.phi_exact(nodes(i));
        }
        double l1; double l2; double linf;
        l1 = l1Norm(sol, sol_exact, grid);
        l2 = l2Norm(sol, sol_exact, grid);
        linf = linfNorm(sol, sol_exact, grid);
        file << N<<","<<h<<","<<l1<<","<<l2<<","<<linf<<"\n";
        H(j) = h;
        l1_errs(j) = l1;
        l2_errs(j) = l2;
        linf_errs(j) = linf;
        j++;    
    }
    double rate_l1 = convergenceRate(H, l1_errs);
    double rate_l2 = convergenceRate(H, l2_errs);
    double rate_linf = convergenceRate(H, linf_errs);
    std::cout << "Convergence rates: L1: " << rate_l1 <<"; L2: " << rate_l2 << "; Linf: " << rate_linf <<"\n";

}
