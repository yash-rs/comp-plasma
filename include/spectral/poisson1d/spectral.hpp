#pragma once

#include <fftw3.h>

namespace spectral{
    class PoissonSpectralSolver{
        double a_; //left domain point
        double b_; //right domain point
        double L_; // lebgth of domain
        double N_; // no. of nodes 0,1,...,N-1

        public:
            PoissonSpectralSolver(double a, double b, double N);

    };
}