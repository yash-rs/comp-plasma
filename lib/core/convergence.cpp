#include "core/convergence.hpp"
#include <Eigen/Dense>
#include <stdexcept>

namespace core {
    double convergenceRate(const Eigen::VectorXd& h, const Eigen::VectorXd& error){
        if(h.size() != error.size() || h.size() < 2){
            throw(std::invalid_argument("Either size of h and error mismatch or length of h is less than 2"));
        }
        if ((h.array() <= 0).any() || (error.array() <= 0).any()) {
            throw(std::invalid_argument("Values of h and error array should be positive"));
        }
        Eigen::VectorXd logh = h.array().log();
        Eigen::VectorXd loge = error.array().log();
        double h_mean = logh.mean();
        double e_mean = loge.mean();
        logh = logh.array() - h_mean;
        loge = loge.array() - e_mean;
        double rate = logh.dot(loge) / logh.dot(logh);
        return rate;
    }
}