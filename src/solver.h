#pragma once

#include "typedefs.h"
#include "system.h"

template<typename T>
struct Solution
{
    Matrix<T>    K;
    Matrix<T>    V;
    Matrix<real> O;

    real theta;

    bool has_eigenvalues;
    bool has_eigenvectors;

    std::vector<T> eigenvalues;
    std::vector<
        Vector<T>
    > eigenvectors;
};

template<typename T>
class Solver
{
public:
    Solver(System* context);
    virtual ~Solver();

    virtual Solution<T> compute(Basis&, real theta=0)=0;
    virtual Solution<T> compute(Basis&, Solution<T>& cache, real theta=0)=0;
    virtual Solution<T> solve(Basis&, bool eigenvectors=0, real theta=0)=0;
    virtual Solution<T> solve(Basis&, Solution<T>& cache,
                              bool eigenvectors=0, real theta=0)=0;
    virtual void        solve(Solution<T>&, bool eigenvectors=0)=0;

protected:
    System* system_context;
};

template<typename T>
class SolverCPU : public Solver<T>
{
public:
    SolverCPU(System* context);
    ~SolverCPU();

    virtual Solution<T> compute(Basis&, real theta=0);
    virtual Solution<T> compute(Basis&, Solution<T>& cache, real theta=0);
    virtual Solution<T> solve(Basis&, bool eigenvectors=0, real theta=0);
    virtual Solution<T> solve(Basis&, Solution<T>& cache,
                              bool eigenvectors=0, real theta=0);
    virtual void        solve(Solution<T>&, bool eigenvectors=0);

protected:
    real overlap(CorrelatedGaussian&, CorrelatedGaussian&);
    T    kinetic(CorrelatedGaussian&, CorrelatedGaussian&, real over);
    T    gaussian_v(real v0, real r0sq, real over, real cij, real theta);
    real c_ij(CorrelatedGaussian&, CorrelatedGaussian&, uint i, uint j);

    void solve_bisection(Solution<T>&, uint max_iterations, real tolorance);
    void solve_full(Solution<T>&, bool eigenvectors);

    void rotate_if_T_complex(T& value, real theta);
};

template<typename T>
Solver<T>::Solver(System* context)
    : system_context(context)
{}

template<typename T>
Solver<T>::~Solver()
{}

