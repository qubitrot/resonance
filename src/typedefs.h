#pragma once

#include <vector>
#include <complex>
#include <eigen3/Eigen/Dense>

typedef double real;
typedef std::complex<real> complex;

constexpr real pi   = 3.14159265358979323846;
constexpr real hbar = 1;

template<typename T>
using Matrix = Eigen::Matrix<T,-1,-1>;

template<typename T>
using Vector = Eigen::Matrix<T,-1,1>;

typedef struct {
    std::vector<short> particle_ids;
    Matrix<real> widths;
    Matrix<real> trans;
    real         norm;
    int          strain;
} CorrelatedGaussian;

typedef std::vector<CorrelatedGaussian> Basis;

//STL specializations
namespace std {
    template<>
    struct hash<pair<short,short>> {
        size_t operator()(pair<short,short> const& ps) const {
            const size_t a ( hash<short>()(ps.first)  );
            const size_t b ( hash<short>()(ps.second) );
            return a^b;
        }
    };

    template<>
    struct equal_to<pair<short,short>> {
        bool operator()(pair<short,short> const& p1,
                        pair<short,short> const& p2) const {
            if ( (p1.first == p2.first  && p1.second == p2.second) ||
                 (p1.first == p2.second && p1.second == p2.first) )
            {
                return true;
            } else {
                return false;
            }
        }
    };
}
