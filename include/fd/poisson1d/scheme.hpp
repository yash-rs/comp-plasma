#pragma once

#include "fd/grid1d.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace fd::poisson1d{
    class Scheme{
        public:
            virtual ~Scheme() = default;
            virtual double centerCoeff(double h) const = 0;
            virtual double neighborCoeff(double h) const = 0;
            virtual void addInteriorRows(std::vector<Eigen::Triplet<double>>&, const Grid1D&, int firstRow, int lastRow) const = 0;
    };

    class SecondOrderScheme : public Scheme{
        public:
            SecondOrderScheme();
            double centerCoeff(double h) const override;
            double neighborCoeff(double h) const override;
            void addInteriorRows(std::vector<Eigen::Triplet<double>>&, const Grid1D&, int firstRow, int lastRow) const override;
    };
}