#include "fd/stencil.hpp"
#include <Eigen/Dense>
#include <cmath>

namespace fd{
    // returns the weights for a general stencil with offsets and given order m
    Eigen::VectorXd stencilWeights(const Eigen::VectorXd &offsets, int order){
        int m = order;
        const int n = offsets.size();
        if(m >= n){
            throw std::invalid_argument("order should be less than number of offsets");
        }
        Eigen::MatrixXd A(n, n);
        Eigen::VectorXd B(n);
        int i, j;
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
                A(i, j) = std::pow(offsets(j), i);
            }
            if (i == m){
                B(i) = std::tgamma(m+1);
            }
            else{
                B(i) = 0;
            }
        }
        Eigen::VectorXd weights = A.colPivHouseholderQr().solve(B);
        return weights;
    }
}