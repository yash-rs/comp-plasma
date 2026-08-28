#pragma once

#include <functional>

namespace core {
    struct ManufacturedSolution{
        std::function<double(double)> rho; //source term which needs to be solved for
        std::function<double(double)> phi_exact; // analytical PDE solution to be tested against
    };

}