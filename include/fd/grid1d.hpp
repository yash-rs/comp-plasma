#pragma once
#include <Eigen/Dense>

namespace fd{
    class Grid1D{
        double a_; //starting point
        double b_; //end point
        double h_; //mesh size
        int N_; // cell count
        Eigen::VectorXd nodes_;
        public:
            Grid1D(double a, double b, int N);
            double a() const {return a_;}
            double b() const {return b_;}
            double h() const {return h_;}
            int N() const {return N_;}
            const Eigen::VectorXd& nodes() const {return nodes_;}

    };
}