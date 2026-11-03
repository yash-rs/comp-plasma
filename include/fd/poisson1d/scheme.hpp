#pragma once

#include "fd/grid1d.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <utility>
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

    class GenScheme{
        std::vector<int> intOffsets_;
        std::vector<int> leftOffsets_;
        std::vector<int> rightOffsets_;
        // raw weights computed internally
        Eigen::VectorXd intWeights_;
        Eigen::VectorXd leftWeights_;
        Eigen::VectorXd rightWeights_;

        static Eigen::VectorXd computeWeights(const std::vector<int>& offsets);

        public:
            GenScheme(const std::vector<int>& intOffsets,
                     const std::vector<int>& leftOffsets,
                     const std::vector<int>& rightOffsets);

            std::pair<std::vector<int>, Eigen::VectorXd> nodeWeights(int i, const Grid1D& grid) const;
            std::pair<std::vector<int>, Eigen::VectorXd> intOffsetsAndWeights() const;
    };
}