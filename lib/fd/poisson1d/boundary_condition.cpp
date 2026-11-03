#include "fd/poisson1d/boundary_condition.hpp"
#include "fd/poisson1d/scheme.hpp"
#include "fd/grid1d.hpp"
#include "fd/stencil.hpp"
#include <complex>
#include <utility>
#include <vector>

namespace fd::poisson1d{
    DirichletBC::DirichletBC(double alpha, double beta) : alpha_(alpha), beta_(beta) {}

    Eigen::SparseMatrix<double> DirichletBC::assembleA(const Grid1D& grid, const Scheme& scheme) const{
        const int n = grid.N() - 1; //interior points for N-1 x N-1 matrix
        const double h = grid.h(); //mesh size
        Eigen::SparseMatrix<double> A(n, n);
        std::vector<Eigen::Triplet<double>> triplets;
        scheme.addInteriorRows(triplets, grid, 1, n-2);
        //row 0
        triplets.emplace_back(0, 0, scheme.centerCoeff(h));
        triplets.emplace_back(0, 1, scheme.neighborCoeff(h));
        //row n-1
        triplets.emplace_back(n-1, n-1, scheme.centerCoeff(h));
        triplets.emplace_back(n-1, n-2, scheme.neighborCoeff(h));
        A.setFromTriplets(triplets.begin(), triplets.end());
        A.makeCompressed();

        return A;
    }

    Eigen::VectorXd DirichletBC::assembleRHS(const Grid1D& grid, const Scheme& scheme, std::function<double(double)> rho) const{
        const int n = grid.N() - 1;
        const double h = grid.h();
        const Eigen::VectorXd nodes = grid.nodes();
        Eigen::VectorXd rhs(n);
        for (int i = 0; i<n; i++){
            rhs(i) = rho(nodes(i+1));
            if(i == 0) rhs(i) += alpha_*-scheme.neighborCoeff(h);
            if(i == n-1) rhs(i) += beta_*-scheme.neighborCoeff(h);
        }

        return rhs;
    }

    Eigen::VectorXd DirichletBC::reconstructSolution(const Grid1D& grid, const Eigen::VectorXd& reduced_sol) const{
        const int N = grid.N();
        Eigen::VectorXd sol(N+1);
        sol(0) = alpha_;
        sol(N) = beta_;
        for (int i=1;i<N;i++){
            sol(i) = reduced_sol(i-1);
        }

        return sol;
    }

    MixedDirichletNeumannBC::MixedDirichletNeumannBC(double alpha, double beta, const std::vector<int>& rightBoundaryOffsets) : alpha_(alpha), beta_(beta), rbo_(rightBoundaryOffsets){}
    Eigen::SparseMatrix<double> MixedDirichletNeumannBC::assembleA(const Grid1D& grid, const GenScheme& scheme) const{
        const int n = grid.N(); // number of unknowns
        const double h = grid.h(); 
        Eigen::SparseMatrix<double> A(n, n);
        std::vector<Eigen::Triplet<double>> triplets;
        Eigen::VectorXd weights;
        std::vector<int> offsets;
        int k;
        std::pair<std::vector<int>, Eigen::VectorXd> res;
        for (int i=0;i<n-1;i++){
            res = scheme.nodeWeights(i+1, grid);
            offsets = res.first;
            weights = res.second;
            k = offsets.size();
            for (int j=0;j<k;j++){
                if(i+1+offsets[j] == 0){
                    continue;
                }
                else{
                    triplets.emplace_back(i, i+offsets[j], weights[j]); // user needs to make sure that offsets dont go out of the grid
                }
            }
        }
        //last row
        int i = n-1;
        k = rbo_.size();
        Eigen::VectorXd o(k);
        for (int i=0;i<k;i++){
            o(i) = rbo_[i];
        }
        weights = stencilWeights(o, 1);
        for(int j=0;j<k;j++){
            triplets.emplace_back(i, i+rbo_[j], (1/h)*weights[j]);
        }
        A.setFromTriplets(triplets.begin(), triplets.end());
        A.makeCompressed();
        return A;
    }

    Eigen::VectorXd MixedDirichletNeumannBC::assembleRHS(const Grid1D& grid, const GenScheme& scheme, std::function<double(double)> rho) const{
        const int n = grid.N();
        Eigen::VectorXd rhs(n);
        const Eigen::VectorXd nodes = grid.nodes();
        Eigen::VectorXd weights;
        std::vector<int> offsets;
        int k;
        std::pair<std::vector<int>, Eigen::VectorXd> res;
        for (int i=0;i<n;i++){
            if (i<=n-2){
                rhs(i) = rho(nodes(i+1));
                res = scheme.nodeWeights(i+1, grid);
                offsets = res.first;
                weights = res.second;
                k = offsets.size();
                for (int j=0;j<k;j++){
                    if(i+1+offsets[j] == 0){
                        rhs(i) -= weights[j]*alpha_;
                    }
                }
            }  
            if (i == n-1){
                rhs(i) = beta_;
            }

        }
        return rhs;
    }

    Eigen::VectorXd MixedDirichletNeumannBC::reconstructSolution(const Grid1D& grid, const Eigen::VectorXd& reduced_sol) const{
        const int N = grid.N();
        Eigen::VectorXd sol(N+1);
        sol(0) = alpha_;
        for (int i=1;i<=N;i++){
            sol(i) = reduced_sol(i-1);
        }
        return sol;
    }

    // ================ Periodic BC =================
    PeriodicBC::PeriodicBC(){}
    Eigen::VectorXd PeriodicBC::assembleRHS(const Grid1D& grid,  std::function<double(double)> rho) const{
        const int n = grid.N();
        Eigen::VectorXd rhs(n);
        const Eigen::VectorXd nodes = grid.nodes();
        for(int i=0;i<n;i++){
            rhs(i) = rho(nodes(i));
        }
        return rhs;
    }
    
    Eigen::VectorXcd PeriodicBC::computeEigs(const Grid1D& grid, const GenScheme& scheme) const{
        const double h = grid.h();
        int n = grid.N();
        double scale = -1 / (h*h);
        std::pair<std::vector<int>, Eigen::VectorXd> res = scheme.intOffsetsAndWeights();
        std::vector<int> o = res.first; // we assume that offsets are symmetric and include 0
        Eigen::VectorXd w = res.second;
        Eigen::VectorXd circulants = Eigen::VectorXd::Zero(n);
        int k = o.size();
        for (int i=0;i<k;i++){
            if(o[i]>=0){
                circulants(o[i]) = w[i];
            }
            else{
                circulants(n+o[i]) = w[i];
            }
            
        }
        Eigen::VectorXcd eigs = Eigen::VectorXcd::Zero(n/2+1);
        for(int m=0;m<=n/2;m++){
            for(int i=0;i<k;i++){
                double theta = 2.0*M_PI*m*o[i]/n;
                eigs(m) += w[i]*std::complex<double>(std::cos(theta), std::sin(theta));
            }
        }
        return eigs;
    }
}