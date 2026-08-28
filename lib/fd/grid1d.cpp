#include "fd/grid1d.hpp"
#include <stdexcept>

namespace fd{
    Grid1D::Grid1D(double a, double b, int N)
    : a_(a), b_(b), h_((b-a)/N), N_(N)
    {
        if (b <= a){
            throw std::invalid_argument("Grid1D: requires b > a");
        }
        if (N < 1){
            throw std::invalid_argument("Grid1D: requires N >= 1");
        }
        nodes_.resize(N_+1);
        for (int i=0; i<=N; i++){
            nodes_(i) = a_+i*h_;
        }
    }
}